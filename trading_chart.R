################# PRICE CHART WEB APP INPIRED BY TRADINGVIEW ###################

library(shiny)
library(htmlwidgets)
library(plotly)
library(data.table)
library(binancer)
library(TTR)
library(plyr)

source('https://raw.githubusercontent.com/Dailer/Astro/master/instruments.R')
instruments=sort(unique(c(bitso_symb0, etoro_symb0, bin_symb0)))

bn_names=c('COMPUSDT','PEPEUSDT','GRTUSDT','IMXUSDT','RONINUSDT','AXLUSDT',
           'SUIUSDT','LUNA2USDT','SGBUSDT','STXUSDT','UNIUSDT','TIAUSDT',
           'ALTUSDT','GMTUSDT','GMXUSDT','HIFIUSDT','LINAUSDT',
           'MASKUSDT','NFPUSDT','PIXELUSDT','PORTALUSDT','PROSUSDT',
           'SAGAUSDT','SUPERUSDT','VIDTUSDT','ZKUSDT','ZROUSDT','ARBUSDT')
yb_names=c('COMP5692USDT','PEPE24478USDT','GRT6719USDT','IMX10603USDT',
           'RON14101USDT','AXL17799USDT','SUI20947USDT','LUNA20314USDT',
           'SGB12186USDT','STX4847USDT','UNI7083USDT','TIA22861USDT',
           'ALT29073USDT','GMT18069USDT','GMX11857USDT','HIFI23037USDT',
           'LINA7102USDT','MASK8536USDT','NFP28778USDT','PIXEL29335USDT',
           'PORTAL29555USDT','PROS8255USDT','SAGA30372USDT','SUPER8290USDT',
           'VIDT22710USDT','ZK24091USDT','ZRO26997USDT','ARB11841USDT')

##### ui -----------------------------------------------------------------------
ui <- fluidPage(
  
  fluidRow(column(4, selectInput("source", "Source",
                              c('Binance' = "binance", 'Yahoo Finance' = "yahoo"), 
                              "binance", width='60%')),
           column(4, selectInput("symb", "Symbol", instruments, 'BTCUSDT',
                                 width='50%')),
           column(4, numericInput('ncand', 'Candles', 100, 60, 940, width='30%'))),
  #column(3, downloadButton('download', 'Download'))
  radioButtons('temp', NULL, c('5m'='5m','15m'='15m','30m'='30m',
                               '1h'='1h','2h'='2h','3h'='3h','4h'='4h',
                               '5h'='5h','6h'='6h','8h'='8h','10h'='10h',
                               '12h'='12h','15h'='15h','18h'='18h','20h'='20h',
                               '1d'='1d','2d'='2d','3d'='3d','5d'='5d','1w'='1w',
                               '1M'='1M'), 
               selected='1h', inline=T),
  hr(),
  plotlyOutput("graph")
)

##### server -------------------------------------------------------------------
server <- function(input, output, session){
  
  # Get kline/candlestick data from Binance, for other timeframes
  binance_klines2=function(symbol, interval=c('10m','20m','25m','45m','3h','5h',
                                              '7h','9h','10h','15h','18h','20h',
                                              '2h','5d'), 
                           limit){
    btf=switch(interval, '10m'='5m', '20m'='5m', '25m'='5m', '45m'='15m', 
               '3h'='1h', '5h'='1h', '7h'='1h', '9h'='1h', '10h'='2h', 
               '15h'='1h','18h'='6h','20h'='4h','2d'='1d','4d'='1d','5d'='1d')
    n1=as.integer(gsub('[[:alpha:]]','',interval))
    n2=as.integer(gsub('[[:alpha:]]','',btf))
    n=n1/n2
    if(limit*n>1000){
      warning('Number of candles must be lesser than ',ceiling(1000/n))
      limit=ceiling(1000/n)-1
    } 
    k=binance_klines(symbol = symbol, interval = btf, limit = limit*n)
    k$id=rep(1:limit, each=n)
    sp=split(k, by='id')
    fmt='%Y-%m-%d %H:%M:%S'
    open_time=unlist(lapply(sp, function(x) strftime(x$open_time[1], format = fmt)))
    close_time=unlist(lapply(sp, function(x) strftime(x$close_time[n], format = fmt)))
    open_time=strptime(open_time, format = fmt)
    close_time=strptime(close_time, format = fmt)
    op=unlist(lapply(sp, function(x) x$open[1]))
    cl=unlist(lapply(sp, function(x) x$close[n]))
    hi=unlist(lapply(sp, function(x) max(x$high)))
    lo=unlist(lapply(sp, function(x) min(x$low)))
    vol=unlist(lapply(sp, function(x) sum(x$volume)))
    options(warn=-1)
    data=data.table(open_time=open_time,open=op,high=hi,low=lo,close=cl,
                    volume=vol,close_time=close_time,symbol=k$symbol[1])
    options(warn=0)
    return(data)
  }
  
  # Retrieve financial data from yahoo finances API webpage
  get_quote_data=function(symbol, interval='1m', range='1d'){
    # interval can be 1m, 2m, 5m, 10m, 15m, 30m, 1h, 2h, 3h, 
    #                 4h, 5h, 6h, 8h, 10h, 12h, 15h, 18h, 20h, 
    #                 1d, 3d, 4d, 5d
    #1m, 2m, 5m, 15m, 30m, 60m, 90m, 1h, 1d, 5d, 1wk, 1mo, 3mo
    require(rvest)
    require(jsonlite)
    require(data.table)
    new_int=c('10m','2h','3h','4h','5h','6h','8h','10h','12h','15h','18h',
              '20h','2d','3d','5d')
    intervs=c('1m','2m','5m','15m','30m','90m','1h','1d','1wk','1mo',
              '3mo', new_int)
    if(!interval %in% intervs) stop('A valid interval must be introduced')
    limit=n=NULL
    # adding new intervals
    if(interval %in% new_int){
      base_int=switch(interval, '10m'='5m', '2h'='1h', '3h'='1h', '4h'='1h', 
                      '5h'='1h', '6h'='1h', '8h'='1h', '10h'='1h', '12h'='1h',
                      '15h'='1h','18h'='1h', '20h'='1h', '3d'='1d', '2d'='1d', 
                      '4d'='1d', '5d'='1d')
      n1=as.integer(gsub('[[:alpha:]]','',interval))
      n2=as.integer(gsub('[[:alpha:]]','',base_int))
      n=n1/n2
      interval=base_int
    }
    # allowing download using a limit candle data
    if(is.numeric(range)){
      limit=as.integer(range)
      if(!is.null(n)) limit=n*limit
      ispl=strsplit(gsub("([a-z]+)",",\\1",interval),",")[[1]]
      range=paste0(as.integer(ispl[1])*limit,ispl[2])
    }
    # scrapping yahoo finances API data
    base_link='https://query1.finance.yahoo.com/v8/finance/chart/'
    link=paste0(base_link,symbol,'?range=',range,'&interval=',interval)
    data=link %>% read_html() %>% html_elements('body') %>% html_text() %>%
      fromJSON()
    # organizing the data in the OHLCV format
    otime=as.POSIXct(data$chart$result$timestamp[[1]], origin="1970-01-01")
    prices=data$chart$result$indicators$quote[[1]]
    res=data.table(open_time=otime, open=prices$open[[1]], high=prices$high[[1]],
                   low=prices$low[[1]], close=prices$close[[1]],
                   volume=prices$volume[[1]], symbol=symbol)
    if(!is.null(limit)) res=tail(res, limit)
    # creating the clandle data for the new intervals
    if(!is.null(n)){
      res$id=rep(1:(limit/n), each=n)
      sp=split(res, by='id')
      otime=sapply(sp, function(x) x$open_time[1])
      otime=as.POSIXct(otime, origin="1970-01-01")
      open=sapply(sp, function(x) x$open[1])
      close=sapply(sp, function(x) x$close[n])
      high=sapply(sp, function(x) max(x$high))
      low=sapply(sp, function(x) min(x$low))
      volume=sapply(sp, function(x) sum(x$volume))
      res=data.table(open_time=otime, open, high, low, close, volume,
                     symbol=symbol)
    }
    res=na.omit(res)
    return(res)
  }
  
  # Get indicators and trading signals for the used strategy
  getsig=function(x){
    require(TTR)
    require(data.table)
    adx=ADX(x[,3:5])[,4]
    adx.norm=(adx-23)/max(abs(adx-23),na.rm=T)
    ema9=EMA(x$close,9)
    ema10=EMA(x$close,10)
    ema24=EMA(x$close,24)
    ema55=EMA(x$close,55)
    fmacd=MACD(x$close, 6, 12, 1)[,2]
    edif=100*(ema9/ema24-1)
    edif.norm=edif/max(abs(edif),na.rm=T)
    sdif=100*(SMA(x$close,9)/SMA(x$close,24)-1)
    sdif.norm=sdif/max(abs(sdif),na.rm=T)
    sqh=c(NA,diff(edif))*sign(sdif)
    idx1=sdif>0 & !is.na(sdif)
    idx2=sdif<=0 & !is.na(sdif)
    scol=rep(NA, nrow(x))
    scol[idx1]=ifelse(sqh[idx1]<0, 2, 1)
    scol[idx2]=ifelse(sqh[idx2]<0, -2, -1)
    adx.sig=ifelse(c(NA,diff(adx-23))<0,T,F)
    sqh.sig=sqh<0
    allsig=adx.sig & sqh.sig
    sig=c(NA,diff(allsig))==1
    allsig=allsig*sign(-sdif)
    sig=sig*sign(-sdif)
    tdif=diff(x$open_time[1:2])
    tf=paste0(as.integer(tdif),toupper(substring(attr(tdif, 'units'),1,1)))
    dt=data.table(x[,1:5],symbol=x$symbol,tf=tf,ADX=adx,ADXn=adx.norm,
                  EMA10=ema10,EMA55=ema55,EMAdif=edif.norm,SMAdif=sdif.norm,
                  class=scol,signal=sig,fmacd=fmacd)
    return(dt)
  }
  
  symbol <- reactive({input$symb})
  src <- reactive({input$source})
  candles <- reactive({input$ncand})
  timeframe <- reactive({input$temp})
  
  ##### get API data -----------------------------------------------------------
  data <- reactive({
    tf=timeframe()
    coin=symbol()
    limit=candles()+55
    if(src()=='binance'){
      if(!coin %in% bin_symb0) coin='XXXUSDT'
      # new implemented timeframes
      ntfs=c('10m','20m','25m','45m','3h','5h','7h','9h','10h','15h','18h',
             '20h','2d','4d','5d')
      if(!tf %in% ntfs){
        kl_data=binance_klines(coin, interval=tf, limit=limit)
      }else{
        kl_data=binance_klines2(coin, interval=tf, limit=limit)
      }
    }else if(src()=='yahoo'){
      if(tf=='1w') tf='1wk'
      if(tf=='1M') tf='1mo'
      if(coin %in% bn_names){
        coin=yb_names[match(coin,bn_names)]
      }
      sym=gsub('USDT', '-USD', coin)
      kl_data=get_quote_data(sym, interval=tf, range=limit)
    }
    sig=getsig(kl_data)
    sig$symbol=symbol()
    sig$source=toupper(src())
    
    tail(sig, candles())
  })
  
  ##### price chart ------------------------------------------------------------
  output$graph <- renderPlotly({
    x=data()
    icol=list(line=list(color='#26A69A'), fillcolor='#26A69A')
    dcol=list(line=list(color='#EF5350'), fillcolor='#EF5350')
    tit=paste(x$symbol[1], tolower(x$tf[1]), x$source[1], sep=' · ')
    # if(length(unique(x$class))<4)
    #   message('Must increase the number of candles to display chart properly')
    x$bcol=plyr::mapvalues(x$class, c(1,2,-1,-2), 
                           c('#B2DF8A','#33A02B','#FC9A99','#E3201D'))
    x$SMAdif=x$SMAdif/max(abs(x$SMAdif),na.rm=T)
    x$ADXn=x$ADXn/max(abs(x$ADXn),na.rm=T)
    x$fmacdn=x$fmacd/max(abs(x$fmacd),na.rm=T)
    
    fig1=plot_ly(data=x, x=~open_time, type="candlestick", open=~open, 
                 close=~close, high=~high, low=~low, increasing=icol, 
                 decreasing=dcol, name='') %>% 
      add_lines(x=~open_time, y=~EMA55, line=list(color='#DC2663', width=2), 
                name='EMA 55', hoverinfo="none", inherit=F, showlegend=F,
                hovertext='') %>%
      layout(yaxis=list(title="USDT", fixedrange=F, spikemode='across', 
                        range='auto', side='right', spikecolor='black',
                        spikedash='dot', spikethickness=1, spikesnap='cursor'), 
             xaxis=list(rangeslider=list(visible=F), spikemode='across',
                        spikecolor='black', spikedash='dot', spikethickness=1,
                        spikesnap='cursor'), 
             title=tit, font=list(size=10))
    #Plotting the indicators
    fig2=plot_ly(data=x, x=~open_time, y=~SMAdif, type='bar', name='SQZMOM',
                 color=~bcol, colors=c('#59AC59','#59FF59','#AC5959','#FF5959'),
                 marker=list(line=list(width=1.5)), showlegend=F) %>%
      add_trace(x=~open_time, y=~ADXn, color=I('black'), type='scatter',
                mode='lines+markers', showlegend=F, name='ADX',
                marker=list(color='black', opacity=0)) %>%
      add_trace(x=~open_time, y=~fmacdn, color=I('blue'), type='scatter',
                mode='lines+markers', showlegend=F, name='Fast MACD',
                marker=list(color='black', opacity=0)) %>%
      layout(xaxis=list(title="Date"), yaxis=list(title="", side="right"))
    subplot(fig1, fig2, heights=c(0.7,0.3), nrows=2, shareX=T, titleY=T) 
  })
  
}

shinyApp(ui, server, options = list(launch.browser = TRUE)) #options = list(launch.browser = TRUE)