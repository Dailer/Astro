
flightpromo=function(){
  require(htm2txt, quietly = T)
  options(warn = -1)
  rl=readLines('https://www.melhoresdestinos.com.br/?s=Cuba+passagens')
  options(warn = 0)
  w=which(grepl('hora|horas| dia| dias|semana|semanas', rl))[1]
  if(length(w)==1){
    txt=gsub('há', 'hace', htm2txt(rl[[w]]))
    txt=paste('Hay ofertas de pasajes para Cuba',txt,'!')
  }else{
    txt='No hay ofertas recientes :('
  }
  message(txt)
}

flightpromo()
