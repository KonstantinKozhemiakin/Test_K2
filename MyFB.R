library(?Rfacebook)
library(httr)

oauth_endpoints("facebook")
AppID <- ""
AppSecret <- ""
app <- oauth_app('facebook', AppID, AppSecret)
fb_token <- oauth2.0_token(oauth_endpoints("facebook"), app,
                           scope = "read_insights",
                           type  = "application/x-www-form-urlencoded", 
                           cache = FALSE)
# This saves the token for future use without having to login everytime you want to use the package
save(fb_token, file = "D://Êîñòè/R/fb_token/fb_token.R")

getFriends(token = fb_token, simplify = FALSE)
MyhomeTL = GET("https://api.twitter.com/1.1/statuses/home_timeline.json", sig)
