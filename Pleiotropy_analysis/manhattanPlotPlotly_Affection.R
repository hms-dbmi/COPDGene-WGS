library( plotly )

### Create df of Phenotype, adjusted pvalue and odds ratio
phewas.df = phewasOutput.df[ ,c("Phenotype", "adjPvalue", "OR")]
phewas.df$Phenotype = as.character(phewas.df$Phenotype)
phewas.df$adjPvalue = as.numeric(phewas.df$adjPvalue)
phewas.df$OR        = as.numeric(phewas.df$OR)

# Calculate log of adjusted Pvalue
phewas.df$logP = -log10(phewas.df$adjPvalue)
# define Position in plot of each pheno
phewas.df$Pos = c(1:nrow(phewas.df))

phewas.df$dotSize = 1
for(i in 1:nrow(phewas.df) ){
  if(phewas.df$logP[i] > 10 ){
    if(phewas.df$OR[i] > 0   & phewas.df$OR[i] <= 0.5){
      phewas.df$dotSize[i] = 2
    }
    if(phewas.df$OR[i] > 0.5 & phewas.df$OR[i] < 1){
      phewas.df$dotSize[i] = 1.5
    }
    if(phewas.df$OR[i] > 1   & phewas.df$OR[i] < 2){
      phewas.df$dotSize[i] = 1.5
    }
    if(phewas.df$OR[i] >= 2  & phewas.df$OR[i] < 3){
      phewas.df$dotSize[i] = 2
    }
    if(phewas.df$OR[i] >= 3  & phewas.df$OR[i] < 4){
      phewas.df$dotSize[i] = 2.5
    }
    if(phewas.df$OR[i] >= 4  & phewas.df$OR[i] < 5){
      phewas.df$dotSize[i] = 3
    }
    if(phewas.df$OR[i] >= 5 ){
      phewas.df$dotSize[i] = 3.5
    }
  }else{
    phewas.df$dotSize[i] = 1 
  }
}

# Adjust dot shape according to being risk (upwards) or predictive (downwards) 
phewas.df$label = ""
phewas.df$dotShape = 16
for(i in 1:nrow(phewas.df)){
  if(as.numeric(phewas.df$logP)[i] > 15){
    phewas.df$label[i] = paste0(phewas.df$Phenotype[i], "(", round(phewas.df$OR[i], 1)  ,")")
  }
  if(phewas.df$OR[i] > 1 & as.numeric(phewas.df$logP)[i] > 3){
    phewas.df$dotShape[i] = 17
  }
  if(phewas.df$OR[i] < 1 & as.numeric(phewas.df$logP)[i] > 3){
    phewas.df$dotShape[i] = 6
  }
}

phewas.df$labelPlotly = paste0(as.character(phewas.df$Phenotype), "(", round(phewas.df$OR, 3), ")")
infun = function(x){
  (max(x) + min(x)) / 2
}

f = list(
  family = "Courier New, monospace",
  size   = 18,
  color  = "#7f7f7f"
)
x = list(
  title     = "Phenotype",
  titlefont = f, 
  tickangle = -40, 
  ticktext  = phewas.df$Phenotype
)

y = list(
  title     = "-log10(p)",
  titlefont = f, 
  ticktext  = phewas.df$Phenotype
)


a = list(
  x = phewas.df$Pos,
  y = phewas.df$logP,
  text = phewas.df$label,
  xref = "x",
  yref = "y",
  showarrow = FALSE,
  arrowhead = 7, 
  tickangle = 40
)


hline = function(y = 0, color = "red") {
  list(
    type = "line", 
    x0   = 0, 
    x1   = 1, 
    xref = "paper",
    y0   = y, 
    y1   = y, 
    line = list(color = color)
  )
}


plot_ly(data = phewas.df,  x = phewas.df$Pos, y = phewas.df$logP, hoverinfo = "text", text = phewas.df$labelPlotly)  %>% 
  layout(title = "PheWAS analysis COPD",
         xaxis = x, 
         yaxis = y, 
         showlegend  = FALSE, 
         #shapes      = hline(15), 
         annotations = a, 
         margin      = list(l = 50, r = 50, b = 150, t = 50, pad = 4))

# subplot(
#   add_markers(p, size = ~as.numeric(phewas.df$OR), symbol = ~as.numeric(phewas.df$dotShape), name = "default"))


