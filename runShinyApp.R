message('library paths:\n', paste('... ', .libPaths(), sep='', collapse='\n'))
#the googlechormeportable app path
chrome.portable = file.path(getwd(),
                            'GoogleChromePortable/App/Chrome-bin/chrome.exe')
# a function that leverage google browser to run the shiny app.
launch.browser = function(appUrl, browser.path=chrome.portable) {
  message('Browser path: ', browser.path)
  shell(sprintf('"%s" --app=%s', browser.path, appUrl))
}
shiny::runApp('./shiny/', launch.browser=launch.browser)