
const macros = {
    "\\R": "\\textsf{R}",
    "\\code": "\\texttt"
};

function processMathHTML()
{
    var l = document.getElementsByClassName('reqn');
    for (let e of l) {
	katex.render(e.textContent, e,
		     {
			 throwOnError: false,
			 macros
		     });
    }
    return;
}

