#ifndef RATICATE_UNKNOWNMATRIX_HPP
#define RATICATE_UNKNOWNMATRIX_HPP

#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "SimpleMatrix.hpp"
#include "SparseArraySeed.hpp"
#include <vector>
#include <memory>
#include <string>
#include <stdexcept>

#ifdef RATICATE_PARALLELIZE_UNKNOWN
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#endif

namespace raticate {

template<typename Data, typename Index>
struct UnknownMatrixCore {
    UnknownMatrixCore(Rcpp::RObject seed) :
        original_seed(seed),
        delayed_env(Rcpp::Environment::namespace_env("DelayedArray")),
        dense_extractor(delayed_env["extract_array"]),
        sparse_extractor(delayed_env["extract_sparse_array"])
    {
        // We assume the constructor only occurs on the main thread, so we
        // won't bother locking things up. I'm also not sure that the
        // operations in the initialization list are thread-safe.

        {
            auto base = Rcpp::Environment::base_env();
            Rcpp::Function fun = base["dim"];
            Rcpp::RObject output = fun(seed);
            if (output.sexp_type() != INTSXP) {
                auto ctype = get_class_name(original_seed);
                throw std::runtime_error("'dim(<" + ctype + ">)' should return an integer vector");
            }
            Rcpp::IntegerVector dims(output);
            if (dims.size() != 2 || dims[0] < 0 || dims[1] < 0) {
                auto ctype = get_class_name(original_seed);
                throw std::runtime_error("'dim(<" + ctype + ">)' should contain two non-negative integers");
            }
            nrow = dims[0];
            ncol = dims[1];
        }

        {
            Rcpp::Function fun = delayed_env["is_sparse"];
            Rcpp::LogicalVector is_sparse = fun(seed);
            if (is_sparse.size() != 1) {
                auto ctype = get_class_name(original_seed);
                throw std::runtime_error("'is_sparse(<" + ctype + ">)' should return a logical vector of length 1");
            }
            sparse = (is_sparse[0] != 0);
        }

        {
            Rcpp::Function fun = delayed_env["chunkdim"];
            Rcpp::RObject output = fun(seed);
            needs_chunks = !output.isNULL();
            if (needs_chunks) {
                Rcpp::IntegerVector chunks(output);
                if (chunks.size() != 2 || chunks[0] < 0 || chunks[1] < 0) {
                    auto ctype = get_class_name(original_seed);
                    throw std::runtime_error("'chunkdim(<" + ctype + ">)' should return a vector containing two non-negative integers");
                }
                chunk_nrow = chunks[0];
                chunk_ncol = chunks[1];
            }
        }

        {
            Rcpp::Function fun = delayed_env["colAutoGrid"];
            Rcpp::RObject output = fun(seed);
            Rcpp::IntegerVector spacing = output.slot("spacings");
            if (spacing.size() != 2 || spacing[1] < 0) {
                auto ctype = get_class_name(original_seed);
                throw std::runtime_error("'spacings' slot of 'colAutoGrid(<" + ctype + ">)' should contain two non-negative integers");
            }
            block_ncol = spacing[1];
        }

        {
            Rcpp::Function fun = delayed_env["rowAutoGrid"];
            Rcpp::RObject output = fun(seed);
            Rcpp::IntegerVector spacing = output.slot("spacings");
            if (spacing.size() != 2 || spacing[0] < 0) {
                auto ctype = get_class_name(original_seed);
                throw std::runtime_error("'spacings' slot of 'rowAutoGrid(<" + ctype + ">)' should contain two non-negative integers");
            }
            block_nrow = spacing[0];
        }
    }

public:
    struct UnknownWorkspace : public tatami::Workspace {
        UnknownWorkspace(bool r = true) : byrow(r) {}
        bool byrow;

        size_t primary_block_start, primary_block_end;
        size_t secondary_chunk_start, secondary_chunk_end;

        std::shared_ptr<tatami::Matrix<Data, Index> > buffer = nullptr;
        std::shared_ptr<tatami::Workspace> bufwork = nullptr;

        Rcpp::RObject contents;
    };

private:
    static Rcpp::RObject create_index_vector(size_t first, size_t last, size_t max) {
        if (first != 0 || last != max) {
            Rcpp::IntegerVector alt(last - first);
            std::iota(alt.begin(), alt.end(), first + 1); // 1-based.
            return alt;
        } else {
            return R_NilValue;
        }
    }

    template<bool byrow>
    Rcpp::List create_quick_indices(size_t i, size_t first, size_t last) const {
        Rcpp::List indices(2);
        indices[(byrow ? 0 : 1)] = Rcpp::IntegerVector::create(i + 1);
        indices[(byrow ? 1 : 0)] = create_index_vector(first, last, (byrow ? ncol : nrow));
        return indices;
    }

    static std::pair<size_t, size_t> round_indices(size_t first, size_t last, size_t interval, size_t max) {
        size_t new_first = (first / interval) * interval;
        size_t new_last = std::min(
            max, 
            (last ? 
                ((last - 1) / interval + 1) * interval // i.e., ceil(last/interval) * interval.
                : 0 
            )
        );
        return std::make_pair(new_first, new_last);
    }

    template<bool byrow>
    Rcpp::List create_rounded_indices(size_t i, size_t first, size_t last, UnknownWorkspace* work) const {
        Rcpp::List indices(2);
        if constexpr(byrow) {
            auto row_rounded = round_indices(i, i + 1, block_nrow, nrow);
            indices[0] = create_index_vector(row_rounded.first, row_rounded.second, nrow);
            work->primary_block_start = row_rounded.first;
            work->primary_block_end = row_rounded.second;

            auto col_rounded = (needs_chunks ? round_indices(first, last, chunk_ncol, ncol) : std::make_pair(first, last));
            indices[1] = create_index_vector(col_rounded.first, col_rounded.second, ncol);
            work->secondary_chunk_start = col_rounded.first;
            work->secondary_chunk_end = col_rounded.second;

        } else {
            auto row_rounded = (needs_chunks ? round_indices(first, last, chunk_nrow, nrow) : std::make_pair(first, last));
            indices[0] = create_index_vector(row_rounded.first, row_rounded.second, nrow);
            work->secondary_chunk_start = row_rounded.first;
            work->secondary_chunk_end = row_rounded.second;

            auto col_rounded = round_indices(i, i + 1, block_ncol, ncol);
            indices[1] = create_index_vector(col_rounded.first, col_rounded.second, ncol);
            work->primary_block_start = col_rounded.first;
            work->primary_block_end = col_rounded.second;
        }
        return indices;
    }

public:
    static bool needs_reset(size_t i, size_t first, size_t last, const UnknownWorkspace* work) {
        bool reset = true;
        if (work->buffer != nullptr) {
            if (i >= work->primary_block_start && i < work->primary_block_end) {
                if (first >= work->secondary_chunk_start && last <= work->secondary_chunk_end) {
                    reset = false;
                }
            }
        }
        return reset;
    }

    template<class Object>
    void check_quick_dense_dims(const Object& obj, size_t first, size_t last) const {
        if (obj.size() != last - first) {
            auto ctype = get_class_name(original_seed);
            throw std::runtime_error("'extract_array(<" + ctype + ">)' returns incorrect dimensions");
        }
    }

    template<bool byrow>
    void quick_dense_extractor_raw(size_t i, Data* buffer, size_t first, size_t last) const {
        auto indices = create_quick_indices<byrow>(i, first, last);
        Rcpp::RObject val0 = dense_extractor(original_seed, indices);

        if (val0.sexp_type() == LGLSXP) {
            Rcpp::LogicalVector val(val0);
            check_quick_dense_dims(val, first, last);
            std::copy(val.begin(), val.end(), buffer);

        } else if (val0.sexp_type() == INTSXP) {
            Rcpp::IntegerVector val(val0);
            check_quick_dense_dims(val, first, last);
            std::copy(val.begin(), val.end(), buffer);

        } else {
            Rcpp::NumericVector val(val0);
            check_quick_dense_dims(val, first, last);
            std::copy(val.begin(), val.end(), buffer);
        }
    }

    template<bool byrow, bool sparse>
    void check_buffered_dims(const tatami::Matrix<Data, Index>* parsed, const UnknownWorkspace* work) const {
        size_t parsed_primary = (byrow ? parsed->nrow() : parsed->ncol());
        size_t parsed_secondary = (byrow ? parsed->ncol() : parsed->nrow());
        size_t expected_primary = work->primary_block_end - work->primary_block_start;
        size_t expected_secondary = work->secondary_chunk_end - work->secondary_chunk_start;

        if (parsed_primary != expected_primary || parsed_secondary != expected_secondary) {
            auto ctype = get_class_name(original_seed);
            throw std::runtime_error("'" + 
                (sparse ? std::string("extract_sparse_array") : std::string("extract_array")) + 
                "(<" + ctype + ">)' returns incorrect dimensions");
        }
    }

    template<bool byrow>
    void buffered_dense_extractor_raw(size_t i, size_t first, size_t last, UnknownWorkspace* work) const {
        auto indices = create_rounded_indices<byrow>(i, first, last, work);
        Rcpp::RObject val0 = dense_extractor(original_seed, indices);

        auto parsed = parse_simple_matrix<Data, Index>(val0);
        check_buffered_dims<byrow, false>(parsed.matrix.get(), work);

        work->buffer = parsed.matrix;
        work->contents = parsed.contents;
        work->bufwork = (work->buffer)->new_workspace(byrow);
    }

public:
    template<class Object>
    void check_quick_sparse_dims(const Object& obj, size_t nnzero) const {
        if (obj.size() != nnzero) {
            auto ctype = get_class_name(original_seed);
            throw std::runtime_error("'extract_sparse_array(<" + ctype + ">)' returns 'nzdata' of the wrong length");
        }
    }

    template<bool byrow>
    void quick_sparse_extractor_raw(size_t i, size_t* n, Data* vbuffer, Index* ibuffer, size_t first, size_t last) const {
        auto indices = create_quick_indices<byrow>(i, first, last);
        Rcpp::RObject val0 = sparse_extractor(original_seed, indices);

        auto dims = parse_dims(val0.slot("dim"));
        int NR = dims.first;
        int NC = dims.second;
        int primary = (byrow ? NR : NC);
        int secondary = (byrow ? NC : NR);
        if (primary != 1 || secondary != static_cast<int>(last - first)) {
            auto ctype = get_class_name(original_seed);
            throw std::runtime_error("'extract_sparse_array(<" + ctype + ">)' returns incorrect dimensions");
        }

        {
            Rcpp::IntegerMatrix indices(Rcpp::RObject(val0.slot("nzindex")));
            if (indices.ncol() != 2) {
                auto ctype = get_class_name(original_seed);
                throw std::runtime_error("'extract_sparse_array(<" + ctype + ">)' should return 'nzindex' with two columns"); 
            }

            *n = indices.rows();
            auto prim = indices.column(byrow ? 0 : 1);
            for (auto p : prim) {
                if (p != 1) {
                    auto ctype = get_class_name(original_seed);
                    throw std::runtime_error("'extract_sparse_array(<" + ctype + ">)' should returns out-of-range 'nzindex'");
                }
            }

            auto idx = indices.column(byrow ? 1 : 0);
            auto icopy = ibuffer;
            for (auto ix : idx) {
                if (ix < 1 || ix > secondary) {
                    auto ctype = get_class_name(original_seed);
                    throw std::runtime_error("'extract_sparse_array(<" + ctype + ">)' should returns out-of-range 'nzindex'");
                }
                *icopy = ix + first - 1; // 0-based indices.
                ++icopy;
            }
        }

        Rcpp::RObject data = val0.slot("nzdata");
        if (data.sexp_type() == LGLSXP) {
            Rcpp::LogicalVector val(data);
            check_quick_sparse_dims(val, *n);
            std::copy(val.begin(), val.end(), vbuffer);

        } else if (data.sexp_type() == INTSXP) {
            Rcpp::IntegerVector val(data);
            check_quick_sparse_dims(val, *n);
            std::copy(val.begin(), val.end(), vbuffer);

        } else {
            Rcpp::NumericVector val(data);
            check_quick_sparse_dims(val, *n);
            std::copy(val.begin(), val.end(), vbuffer);
        }

        return;
    }

    template<bool byrow>
    void buffered_sparse_extractor_raw(size_t i, size_t first, size_t last, UnknownWorkspace* work) const {
        auto indices = create_rounded_indices<byrow>(i, first, last, work);
        auto val0 = sparse_extractor(original_seed, indices);

        auto parsed = parse_SparseArraySeed<Data, Index>(val0);
        check_buffered_dims<byrow, true>(parsed.matrix.get(), work);

        work->buffer = parsed.matrix;
        work->contents = parsed.contents;
        work->bufwork = (work->buffer)->new_workspace(byrow);
    }

public:
    size_t nrow, ncol;
    bool sparse;

private:
    bool needs_chunks;
    size_t chunk_nrow, chunk_ncol;
    size_t block_nrow, block_ncol;

    Rcpp::RObject original_seed;
    Rcpp::Environment delayed_env;
    Rcpp::Function dense_extractor, sparse_extractor;
};

#ifdef RATICATE_PARALLELIZE_UNKNOWN

/*****************************
 *** Main thread evaluator ***
 *****************************/

template<typename Data, typename Index>
struct UnknownEvaluator {
    bool sparse;
    bool buffered;
    bool byrow;

    size_t index, first, last;
    size_t* nonzero;
    Data* dbuffer;
    Index* ibuffer;
    typename UnknownMatrixCore<Data, Index>::UnknownWorkspace* work;

    const UnknownMatrixCore<Data, Index>* parent;

    bool parallel = false;
    bool ready = false;
    bool finished = false;
    std::string error;

    bool create_work = false;
    typename UnknownMatrixCore<Data, Index>::UnknownWorkspace** new_work;

public:
    template<bool B>
    void set(size_t i, size_t f, size_t l, const UnknownMatrixCore<Data, Index>* core) {
        byrow = B;
        index = i;
        first = f;
        last = l;
        parent = core;
        ready = true;
        finished = false;
        create_work = false;
    }

    template<bool B>
    void set(size_t i, Data* buffer, size_t f, size_t l, const UnknownMatrixCore<Data, Index>* core) {
        set<B>(i, f, l, core);
        sparse = false;
        buffered = false;
        dbuffer = buffer;
    }

    template<bool B>
    void set(size_t i, Data* buffer, size_t f, size_t l, typename UnknownMatrixCore<Data, Index>::UnknownWorkspace* w, const UnknownMatrixCore<Data, Index>* core) {
        set<B>(i, f, l, core);
        sparse = false;
        buffered = true;
        dbuffer = buffer;
        work = w;
    }

    template<bool B>
    void set(size_t i, size_t* n, Data* data, Index* index, size_t f, size_t l, const UnknownMatrixCore<Data, Index>* core) {
        set<B>(i, f, l, core);
        sparse = true;
        buffered = false;
        nonzero = n;
        dbuffer = data;
        ibuffer = index; 
    }

    template<bool B>
    void set(size_t i, Data* data, Index* index, size_t f, size_t l, typename UnknownMatrixCore<Data, Index>::UnknownWorkspace* w, const UnknownMatrixCore<Data, Index>* core) {
        set<B>(i, f, l, core);
        sparse = true;
        buffered = true;
        dbuffer = data;
        ibuffer = index; 
        work = w;
    }

public:
    void set(typename UnknownMatrixCore<Data, Index>::UnknownWorkspace** nw, bool B) {
        new_work = nw;
        byrow = B;
        create_work = true;
        ready = true;
        finished = false;
    }

public:
    void harvest() {
        if (create_work) {
            *new_work = new typename UnknownMatrixCore<Data, Index>::UnknownWorkspace(byrow);
        } else {
            if (!sparse) {
                if (buffered) {
                    if (byrow) {
                        parent->template buffered_dense_extractor_raw<true>(index, first, last, work);
                    } else {
                        parent->template buffered_dense_extractor_raw<false>(index, first, last, work);
                    }
                } else {
                    if (byrow) {
                        parent->template quick_dense_extractor_raw<true>(index, dbuffer, first, last);
                    } else {
                        parent->template quick_dense_extractor_raw<false>(index, dbuffer, first, last);
                    }
                }
            } else {
                if (buffered) {
                    if (byrow) {
                        parent->template buffered_sparse_extractor_raw<true>(index, first, last, work);
                    } else {
                        parent->template buffered_sparse_extractor_raw<false>(index, first, last, work);
                    }
                } else {
                    if (byrow) {
                        parent->template quick_sparse_extractor_raw<true>(index, nonzero, dbuffer, ibuffer, first, last);
                    } else {
                        parent->template quick_sparse_extractor_raw<false>(index, nonzero, dbuffer, ibuffer, first, last);
                    }
                }
            }
        }
        finished = true;
    }

    void reset() {
        finished = false;
        ready = false;
    }
};

template<typename Data, typename Index>
UnknownEvaluator<Data, Index>& unknown_evaluator() {
    static UnknownEvaluator<Data, Index> e;
    return e;
}

/****************************
 *** Parallel coordinator ***
 ****************************/

struct ParallelCoordinator {
    std::mutex coord_lock;
    std::mutex rcpp_lock;
    std::condition_variable cv;

    template<typename Data, typename Index>
    struct OnMainExit {
        UnknownEvaluator<Data, Index> copy;
        OnMainExit() : copy(unknown_evaluator<Data, Index>()) {}
        ~OnMainExit() {
            auto& ex = unknown_evaluator<Data, Index>();
            ex = copy;
        }
    };

    template<typename Data, typename Index, class Function>
    void run(size_t n, Function f, size_t nthreads) {
        // Acquire the evaluator lock to indicate that we're currently in a single
        // parallel context. This avoids wacky messages from other calls to run().
        std::lock_guard<std::mutex> clk(coord_lock);

        // This restores the state of the unknown evaluator to what it was
        // before entry into this function, to enable nested calls to run().
        auto& ex = unknown_evaluator<Data, Index>();
        OnMainExit<Data, Index> copier;

        // Only parallelize if it's strictly necessary.
        ex.parallel = (n > 1 && nthreads > 1);
        ex.error = "";

        if (ex.parallel) {
            size_t jobs_per_worker = std::ceil(static_cast<double>(n) / nthreads);
            size_t start = 0;
            std::vector<std::thread> jobs;
            std::atomic_size_t ncomplete = 0;
            std::vector<std::string> errors(nthreads);

            for (size_t w = 0; w < nthreads; ++w) {
                size_t end = std::min(n, start + jobs_per_worker);
                if (start >= end) {
                    ncomplete++;
                    continue;
                }

                // Local scope, avoid shenanigans when 'w' increments.
                size_t id = w;

                jobs.emplace_back([&](size_t s, size_t e) -> void {
                    try {
                        f(s, e);
                    } catch (std::exception& x) {
                        // No throw here, we need to make sure we mark the
                        // thread as being completed so that the main loop can quit.
                        errors[id] = x.what();
                    }
                    ncomplete++;
                    cv.notify_all();
                }, start, end);

                start += jobs_per_worker;
            }

            // Handling all requests from the workers.
            while (1) {
                std::unique_lock lk(rcpp_lock);
                cv.wait(lk, [&]{ return (ex.ready && !ex.finished) || ncomplete.load() == nthreads; });
                if (ncomplete.load() == nthreads) {
                    break;
                }

                try {
                    ex.harvest();
                } catch (std::exception& x) {
                    // No throw, we need to make sure we notify the worker of (failed) completion.
                    ex.finished = true;
                    ex.error = x.what();
                } catch (...) {
                    // Sometimes R throws these weird Rcpp::LongjumpException errors.
                    ex.finished = true;
                    ex.error = "failed extraction from the unknown matrix";
                }

                // Unlock before notifying, see https://en.cppreference.com/w/cpp/thread/condition_variable
                lk.unlock();
                cv.notify_all();
            }

            for (auto& job : jobs) {
                job.join();
            }

            for (auto err : errors) {
                if (!err.empty()) {
                    throw std::runtime_error(err);
                }
            }
        } else {
            f(0, n);
        }
    }

    template<typename Data, typename Index, class ParallelFunction, class SerialFunction>
    void lock(ParallelFunction pfun, SerialFunction sfun) {
        auto& ex = unknown_evaluator<Data, Index>();
        if (!ex.parallel) {
            sfun();
            return;
        }

        // Waiting until the main thread executor is free,
        // and then assigning it a task.
        {
            std::unique_lock lk(rcpp_lock);
            cv.wait(lk, [&]{ return !ex.ready; });

            // We can throw here because we're not obliged to initiate
            // discussion with the main thread; so it's not like we're
            // going to be leaving the main thread in a hanging state.
            if (!ex.error.empty()) {
                throw std::runtime_error(ex.error);
            }

            try {
                pfun();
            } catch (std::exception& e) {
                // Task assignment failed, so we make sure to reset to
                // avoid partial assignment (and unblock all other workers).
                ex.reset();
                ex.error = e.what();
                throw;
            }
        }

        // Notifying everyone that there is a task. At this point,
        // ready = true and finished = false, so the waiting workers
        // should not proceed; only the main thread should respond.
        cv.notify_all();

        // Checking that we get the finished result, and then we set
        // ready = false to give another worker thread the chance to acquire the lock.
        {
            std::unique_lock lk(rcpp_lock);
            cv.wait(lk, [&]{ return ex.finished; });

            ex.reset();
            if (!ex.error.empty()) {
                // Throwing after the reset so that other workers are not blocked.
                throw std::runtime_error(ex.error);
            }
        }
    }
};

inline ParallelCoordinator& parallel_coordinator() {
    static ParallelCoordinator c;
    return c;
}

#endif

template<typename Data, typename Index>
class UnknownMatrix : public tatami::Matrix<Data, Index> {
public:
    UnknownMatrix(Rcpp::RObject seed) : core(seed) {}

    size_t nrow() const {
        return core.nrow;
    }

    size_t ncol() const {
        return core.ncol;
    }

    bool sparse() const {
        return core.sparse;
    }

    bool prefer_rows() const {
        // All of the individual extract_array outputs are effectively column-major.
        return false;
    }

public:
    std::shared_ptr<tatami::Workspace> new_workspace(bool row) const { 
        std::shared_ptr<tatami::Workspace> output;

#ifdef RATICATE_PARALLELIZE_UNKNOWN 
        // We default-initialize an Rcpp::RObject, so we lock it just in case.
        typename UnknownMatrixCore<Data, Index>::UnknownWorkspace* tmp;
        auto& par = parallel_coordinator();
        auto& ex = unknown_evaluator<Data, Index>();
        par.template lock<Data, Index>(
            [&]() -> void {
                ex.set(&tmp, row);
            },
            [&]() -> void {
                tmp = new typename UnknownMatrixCore<Data, Index>::UnknownWorkspace(row);
            }
        );
        output.reset(tmp);
#else
        output.reset(new typename UnknownMatrixCore<Data, Index>::UnknownWorkspace(row));
#endif

        return output;
    }

private:
    template<bool byrow>
    void quick_dense_extractor(size_t i, Data* buffer, size_t first, size_t last) const {
#ifndef RATICATE_PARALLELIZE_UNKNOWN 
        core.template quick_dense_extractor_raw<byrow>(i, buffer, first, last);
#else
        auto& ex = unknown_evaluator<Data, Index>();
        auto& par = parallel_coordinator();
        par.template lock<Data, Index>(
            [&]() -> void {
                ex.template set<byrow>(i, buffer, first, last, &core);
            },
            [&]() -> void {
                core.template quick_dense_extractor_raw<byrow>(i, buffer, first, last);
            }
        );
#endif
    }

    template<bool byrow>
    void buffered_dense_extractor(size_t i, Data* buffer, size_t first, size_t last, tatami::Workspace* work0) const {
        auto work = static_cast<typename UnknownMatrixCore<Data, Index>::UnknownWorkspace*>(work0);
        if (work->byrow != byrow) {
            throw std::runtime_error("workspace should have been generated with 'row=" + std::to_string(byrow) + "'");
        }

        if (core.needs_reset(i, first, last, work)) {
#ifndef RATICATE_PARALLELIZE_UNKNOWN 
            core.template buffered_dense_extractor_raw<byrow>(i, first, last, work);
#else
            auto& ex = unknown_evaluator<Data, Index>();
            auto& par = parallel_coordinator();
            par.template lock<Data, Index>(
                [&]() -> void {
                    ex.template set<byrow>(i, buffer, first, last, work, &core);
                },
                [&]() -> void {
                    core.template buffered_dense_extractor_raw<byrow>(i, first, last, work);
                }
            );
#endif
        }

        i -= work->primary_block_start;
        first -= work->secondary_chunk_start;
        last -= work->secondary_chunk_start;
        if constexpr(byrow) {
            (work->buffer)->row_copy(i, buffer, first, last, (work->bufwork).get());
        } else {
            (work->buffer)->column_copy(i, buffer, first, last, (work->bufwork).get());
        }
    }

public:
    const Data* row(size_t r, Data* buffer, size_t first, size_t last, tatami::Workspace* work=nullptr) const {
        if (work == NULL) {
            quick_dense_extractor<true>(r, buffer, first, last);
        } else {
            buffered_dense_extractor<true>(r, buffer, first, last, work);
        }
        return buffer;
    }

    const Data* column(size_t c, Data* buffer, size_t first, size_t last, tatami::Workspace* work=nullptr) const {
        if (work == NULL) {
            quick_dense_extractor<false>(c, buffer, first, last);
        } else {
            buffered_dense_extractor<false>(c, buffer, first, last, work);
        }
        return buffer;
    } 

private:
    template<bool byrow>
    tatami::SparseRange<Data, Index> quick_sparse_extractor(size_t i, Data* vbuffer, Index* ibuffer, size_t first, size_t last, bool sorted) const {
        size_t n = 0;

#ifndef RATICATE_PARALLELIZE_UNKNOWN
        core.template quick_sparse_extractor_raw<byrow>(i, &n, vbuffer, ibuffer, first, last);
#else
        auto& ex = unknown_evaluator<Data, Index>();
        auto& par = parallel_coordinator();
        par.template lock<Data, Index>(
            [&]() -> void {
                ex.template set<byrow>(i, &n, vbuffer, ibuffer, first, last, &core);
            },
            [&]() -> void {
                core.template quick_sparse_extractor_raw<byrow>(i, &n, vbuffer, ibuffer, first, last);
            }
        );
#endif

        if (sorted && !std::is_sorted(ibuffer, ibuffer + n)) {
            // TODO: use an in-place sort?
            std::vector<std::pair<Data, Index> > holding;
            holding.reserve(n);
            for (size_t ix = 0; ix < n; ++ix) {
                holding.emplace_back(vbuffer[ix], ibuffer[ix]);
            }
            std::sort(holding.begin(), holding.end());
            for (size_t ix = 0; ix < n; ++ix) {
                vbuffer[ix] = holding[ix].first;
                ibuffer[ix] = holding[ix].second;
            }
        }

        return tatami::SparseRange<Data, Index>(n, vbuffer, ibuffer);
    }

    template<bool byrow>
    tatami::SparseRange<Data, Index> buffered_sparse_extractor(size_t i, Data* vbuffer, Index* ibuffer, size_t first, size_t last, tatami::Workspace* work0, bool sorted) const {
        auto* work = static_cast<typename UnknownMatrixCore<Data, Index>::UnknownWorkspace*>(work0);
        if (work->byrow != byrow) {
            throw std::runtime_error("workspace should have been generated with 'row=" + std::to_string(byrow) + "'");
        }

        if (core.needs_reset(i, first, last, work)) {
#ifndef RATICATE_PARALLELIZE_UNKNOWN
            core.template buffered_sparse_extractor_raw<byrow>(i, first, last, work);
#else
            auto& ex = unknown_evaluator<Data, Index>();
            auto& par = parallel_coordinator();
            par.template lock<Data, Index>(
                [&]() -> void {
                    ex.template set<byrow>(i, vbuffer, ibuffer, first, last, work, &core);
                },
                [&]() -> void {
                    core.template buffered_sparse_extractor_raw<byrow>(i, first, last, work);
                }
            );
#endif
        }

        i -= work->primary_block_start;
        first -= work->secondary_chunk_start;
        last -= work->secondary_chunk_start;

        tatami::SparseRange<Data, Index> output;
        if constexpr(byrow) {
            output = (work->buffer)->sparse_row_copy(i, vbuffer, ibuffer, first, last, tatami::SPARSE_COPY_BOTH, (work->bufwork).get(), sorted);
        } else {
            output = (work->buffer)->sparse_column_copy(i, vbuffer, ibuffer, first, last, tatami::SPARSE_COPY_BOTH, (work->bufwork).get(), sorted);
        }

        // Need to adjust the indices.
        for (size_t i = 0; i < output.number; ++i) {
            ibuffer[i] += work->secondary_chunk_start;
        }

        return output;
    }

public:
    tatami::SparseRange<Data, Index> sparse_row(size_t r, Data* vbuffer, Index* ibuffer, size_t first, size_t last, tatami::Workspace* work=nullptr, bool sorted=true) const {
        if (core.sparse) {
            if (work == NULL) {
                return quick_sparse_extractor<true>(r, vbuffer, ibuffer, first, last, sorted);
            } else {
                return buffered_sparse_extractor<true>(r, vbuffer, ibuffer, first, last, work, sorted);
            }
        } else {
            return tatami::Matrix<Data, Index>::sparse_row(r, vbuffer, ibuffer, first, last, work, sorted);
        }
    }

    tatami::SparseRange<Data, Index> sparse_column(size_t c, Data* vbuffer, Index* ibuffer, size_t first, size_t last, tatami::Workspace* work=nullptr, bool sorted=true) const {
        if (core.sparse) {
            if (work == NULL) {
                return quick_sparse_extractor<false>(c, vbuffer, ibuffer, first, last, sorted);
            } else {
                return buffered_sparse_extractor<false>(c, vbuffer, ibuffer, first, last, work, sorted);
            }
        } else {
            return tatami::Matrix<Data, Index>::sparse_column(c, vbuffer, ibuffer, first, last, work, sorted);
        }
    }

private:
    UnknownMatrixCore<Data, Index> core;
};

}

#endif
