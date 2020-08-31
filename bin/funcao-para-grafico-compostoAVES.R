# Função para GRÁFICO COMPOSTO (Genérico)
### Criada por Victor Lemes Landeiro
######### Última atualização em 18-07-2008 ###

generico<-function(tabela,gradiente,at,grad,eixoY,eixoX){ ## Começo da função onde são especificados todos os 6 argumentos:
  # “tabela”, “gradiente”, at, grad, eixoY, eixoX. Abaixo seguem alguns comentários e exemplos sobre a função e seus argumentos.
  ####### Argumentos
  ########generico<-function(tabela, gradiente ,at,"grad","eixoY","eixoX")
  # tabela é a tabela de espécies, uma data frame ou matriz com os dados das espécies
  # gradiente arquivo com os valores do gradiente
  # at é usado para alterar a posição do nome das espécies no gráfico, comece com o valor 1 e vá aumentando até os nomes ficarem na posição desejada.
  # grad é o texto para a legenda que deseja colocar no gráfico mostrando o nome do gradiente. Esta legenda não pode ser muito longa. Deve vir entre aspas.
  # eixoY e eixoX são as legendas que deseja colocar nos eixos x e Y. Deve vir entre aspas.
  ##Para colocar o nome das espécies no local desejado é preciso mudar o valor do argumento at
  #### Abaixo segue o código da função.
  tabela <- as.matrix(tabela)
  gradiente <- as.matrix(gradiente)
  media.pond <- colSums(tabela*gradiente[,1])/colSums(tabela)
  sub.orden <- tabela[order(gradiente[,1],decreasing=F),]	# Ordenar parcelas de acordo com o gradiente
  sub.orde <- sub.orden[,order(media.pond,decreasing=T)] # colocar espécies ordenadas pela média ponderada
  dados.pa <- matrix(0,nrow(tabela),ncol(tabela))
  dados.pa[tabela>0]<-1
  ordenado <- sub.orde[,which(colSums(dados.pa)>0)] ## para deletar possíveis colunas vazias (espécie que não ocorreu)
  par(mfrow=c(ncol(ordenado)+1,1), mar=c(0,4,0.2,10),oma=c(3,1,1,6))
  layout(matrix(1:(ncol(ordenado)+1)), heights=c(3,rep(1,ncol(ordenado))))
  plot(sort(gradiente[,1]),axes=F,ylab="",mfg=c(21,1),lwd=10,las=2,lend="butt",frame.plot=F,xaxt="n",type="h",col="black",ylim=c(min(gradiente),max(gradiente)))
  axis(side=2,at=c(0,max(gradiente)),las=2)
  mtext(grad,4,outer=T,font=2,line=-10,padj=-18.5,las=2)
  for(i in 1:ncol(ordenado)){
    barplot(ordenado[,i],bty="l",axisnames=F,axes=FALSE,col="black")#
    axis(side=2,at=max(ordenado[,i]),las=2)
    mtext(colnames(ordenado)[i],3,line=-1.0,adj=0,at=at,cex=.8,font=3)
    }
  mtext(eixoX,1,outer=T,font=2,line=1.2)
  mtext(eixoY,2,font=2,outer=T,line=-2)
  }
  # Este gráfico será construído através de uma matriz de dados de espécies (Presença e ausência ou abundância) que será ordenada de acordo com a média ponderada, calculada a partir dos dados de um gradiente ecológico.
  ############## Até Aqui! ###################################