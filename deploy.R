library(rsconnect)
library(shiny)

rsconnect::setAccountInfo(name = "pcbrendel",
                          token = "DE39BA7F2A6142565F38C634F2F85353",
                          secret = "pKLCILPxLRz/l4jjebJqi2zjY1ecIUFrSlGaNgYD")

# test app
runApp()

# deploy app
deployApp()

# terminate
terminateApp("<your app's name>")