#classification demo w/ hpo terms
library(shiny)
library(shinydashboard)
library(Morpho)
library(rgl)
library(ggplot2)
library(dplyr)
library(geomorph)
library(stats)
library(caret)
library(shapes)
library(sparsediscrim)
library(corpcor)
library(plotly)
library(shinycssloaders)
library(Rvcg)
library(shinyjs)

# save(, file = "~/shiny/Classification_demo/demo_objects.Rdata")
# load("/srv/shiny-server/Classification_demo/demo_objects.Rdata")
load("/Users/jovid/shiny/shinyapps/FB2_HPO_classification/app_startup.Rdata")
load("~/Downloads/mshape.Rdata")

body <- dashboardBody(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  fluidRow(
    box(title = tags$b("About this app"),
        solidHeader = T,
        status = "warning",
        collapsible = T,
        width = 12,
        "This application serves as a demonstration of our syndrome classifier. Choose one of the preloaded faces under the [Pick Parameters] tab. The face will be displayed under the [Scan] tab. The [Syndrome Probabilities] tab will show the top ten syndrome probabilities from the classifier. You can compare the individual's facial shape similarity to the average shape of any syndrome by checking the Make a comparison box under the [Pick parameters] tab."
    ),
    box(title = tags$b("Scan"),
        solidHeader = T,
        status = "warning",
        collapsible = F,
        withSpinner(rglwidgetOutput("LM_face", width = "auto"), type = 6, color = "#757471")
        ),
    box(title = tags$b("Pick parameters"),
        solidHeader = T,
        status = "warning",
        collapsible = F,
        selectInput("file1", "", choices = c("Apert", "Achondroplasia", "Williams", "Treacher Collins", "Noonan", "22q deletion", "Cockayne", "Nager"), selected = "22q deletion", multiple = F),
        splitLayout(cellWidths = c("50%", "50%"),
                    checkboxGroupInput("compare", label = "Morphology", choices = "Make a comparison?"),
                    checkboxInput("dense", label = "Dense landmarks")
                    ),
        splitLayout(cellWidths = c("50%", "50%"),
                    checkboxGroupInput("hpo", label = "Extra phenotype info", choices = "Add HPO terms"),
                    disabled(downloadButton("report", "Send me the results!"))
                    ),
        conditionalPanel(condition = "input.compare == 'Make a comparison?'",
                         selectInput("reference", label = "Syndrome", choices = sort(unique(hdrda.df$synd)), selected = "Non-syndromic"), 
                         checkboxInput("displace", label = "Plot displacement vectors?", value = F),
                         conditionalPanel("input.displace == 1",
                                          sliderInput("transparency", label = "Mesh transparency", min = 0, max = 1, step = .1, value = .5)
                         )
        ),
        conditionalPanel(condition = "input.hpo == 'Add HPO terms'",
                         selectInput("hpo_terms", label = "Terms", choices = as.character(hpo$name[names(hpo$name) %in% unique(phenotype_2022$HPO_ID)]), selected = "Pointed chin")
                        )
      )
    ),
    box(title = tags$b("Syndrome probabilities"),
        status = "warning",
        solidHeader = T,
        collapsible = F,
        withSpinner(plotlyOutput("diagnosis_scree"), type = 6, color = "#757471")
        ),
    box(title = tags$b("HPO-related info"),
        status = "warning",
        solidHeader = T,
        collapsible = T,
        collapsed = T,
        dataTableOutput('hpo_synds')
      )
  )


# We'll save it in a variable `ui` so that we can preview it in the console
dbHeader <- dashboardHeader()
dbHeader$children[[2]]$children <-  tags$a(href='../',
                                                     tags$img(src="uc_logo.png",height='30',width='116'))
ui <- dashboardPage(title = "Classification demo",
  dbHeader,
  dashboardSidebar(disable = T),
  body
)

server <- function(input, output, session){
  
  name.key <- c(Apert = 1419, Achondroplasia = 34, Williams = 1702, `Treacher Collins` = 47, Noonan =699, `22q deletion` = 1357, Cockayne = 2423, Nager = 799)
  # 
  # outVar <- reactive({
  #   hpo.terms <- hpo$name[grep(tolower(hpo$name), pattern = tolower(input$hpo_input))]
  #   hpo.id <- names(hpo.terms)
  #   hpo.terms <- as.character(hpo.terms)
  #   return(list(hpo.id, hpo.terms))
  # })
  # 
  # #slow down reactivity of text input
  # ontology.list <- debounce(outVar, 1200)
  # 
  # output$variables <- renderUI({
  #   selectInput('variables2', 'Associated HPO terms', ontology.list()[[2]], multiple = T)
  # })


  mesh.et.lms <- eventReactive(input$file1, {
  # mesh.et.lms <- reactive({
    # file1 <- name.key[names(name.key) %in% "Nager"]
    file1 <- name.key[names(name.key) %in% input$file1]
    
    file.mesh <- atlas
    synd.mshape <- mean.shape
    file.mesh$vb[-4,] <- t(showPC(as.matrix(hdrda.df[file1, -1]), PC.eigenvectors, synd.mshape)[,,1])

    return(list(file.mesh))
  })

  output$LM_face <- renderRglwidget({
    
    clear3d()
    bg3d(color = rgb(245/255, 245/255, 245/255, .9))
    
    # par3d(userMatrix = matrix(c(.998,-.005,.0613,0,.0021,.999,.045,0,-.061,-.045,.997,0,0,0,0,1),ncol =4,nrow = 4))
    par3d(userMatrix = diag(4))
    aspect3d("iso")
    par3d(zoom = .7)
    
    sample1k <- sample(1:27903, 1000)
    meshrot <- rotmesh.onto(file.mesh, t(file.mesh$vb[-4, sample1k]), t(atlas$vb[-4, sample1k]))
    meshrot <- rotmesh.onto(mesh.et.lms()[[1]], t(mesh.et.lms()[[1]]$vb[-4, sample1k]), t(atlas$vb[-4, sample1k]))

    plot3d(vcgSmooth(meshrot$mesh), col = "lightgrey", alpha = 1, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
    # spheres3d(as.matrix(meshrot$yrot), radius = .75, color = "red")

    if(input$dense > 0) points3d(t(meshrot$mesh$vb)[,-4], col = rgb(1, 0, 0), alpha = .1)

    if(length(input$compare) > 0){

      #find individuals' registered lms
      ind.lms <- matrix(as.numeric(filtered.lms[as.character(filtered.lms[,1]) %in% mesh.et.lms()[[3]],-1:-2]), nrow = 65, byrow = T)

      #rotonto atlas
      scaled.lms <- procOPA(mshape(arrayspecs(filtered.lms[,-1:-2], 65, 3)), as.matrix(atlas.lms))$Bhat
      scaled.mesh <- rotmesh.onto(atlas, refmat = as.matrix(atlas.lms), tarmat = as.matrix(scaled.lms), scale = T, reflection = T)

      #calculate specified sydrome mean
      syndrome.mean <- matrix(as.numeric(colMeans(filtered.lms[filtered.lms$Syndrome == input$reference,-1:-2])), nrow = 65, byrow = T)
      #tps3d
      clear3d()

      # par3d(userMatrix = matrix(c(.998,-.005,.0613,0,.0021,.999,.045,0,-.061,-.045,.997,0,0,0,0,1),ncol =4,nrow = 4))
      par3d(userMatrix = matrix(c(-.017,-.999,-.022,0,-.999,-.016,-.03,0,.03,-.023,.999,0,0,0,0,1),ncol =4,nrow = 4))
      par3d(zoom = .7)


      synd.mesh <- tps3d(scaled.mesh$mesh, scaled.lms, syndrome.mean)
      ind.mesh <- rotmesh.onto(mesh.et.lms()[[1]], as.matrix(mesh.et.lms()[[2]]), ind.lms, scale = T, reflection = T)$mesh

      if(input$displace == F){
        bg3d(color = rgb(245/255, 245/255, 245/255, .9))
        mD.synd <- meshDist(ind.mesh, synd.mesh , plot = F, scaleramp = F, displace = input$displace, alpha = 1)
        a <- render(mD.synd, displace = input$displace, alpha = 1)
      } else if(input$displace){
        bg3d(color = "#1a1a1a")
        mD.synd <- meshDist(ind.mesh, synd.mesh, plot = F, scaleramp = F, displace = input$displace, alpha = input$transparency)
        a <- render(mD.synd, displace = input$displace, alpha = input$transparency)
      }

    }
    aspect3d("iso")
    rglwidget()

  })
  # 
  # #change to enable when a scan is uploaded: DONE
  # observeEvent(is.null(input$file1) == F,
  #              {
  #                disable("report")
  #              })
  # 
  # #update the priors if given hpo terms
  # # new.mod <- eventReactive(input$update_model,{
  new.mod <- reactive({
    hpo.term <- unique(phenotype_2022$HPO_ID[phenotype_2022$HPO_ID == names(hpo$name)[which(as.character(hpo$name) == input$hpo_terms)]])   #just our dataset HPO

    #see all hpo names for current synd: 
    # View(data.frame(hpo$name[names(hpo$name) %in% unique(standardized.freqs$term[standardized.freqs$synd == official.names[i]])]))
    # View(data.frame(hpo$name[names(hpo$name) %in% unique(hpo.pos[hpo.pos$V3 == official.names[i],5])]))
    #see current hpo name: hpo$name[names(hpo$name) == hpo.term]
    in.hpo <- rep(NA, length(official.names))
    for(k in 1:length(in.hpo)) in.hpo[k] <- length(grep(official.names[k], x = phenotype_2022$DiseaseName[phenotype_2022$HPO_ID == hpo.term])) > 0
    
    in.hpo[official.names == "Non-syndromic"] <- TRUE
    #debug: levels(hdrda.df$synd)[in.hpo]
    
    #priors are adjusted to uniformly sharing 10% for each syndrome not in the selected HPO, the rest of the weight is uniform with the remaining syndromes in the HPO list
    updated.priors <- rep(NA, length(official.names))
    names(updated.priors) <- official.names
    
    updated.priors[in.hpo == F] <- 0.1 / length(official.names[in.hpo == F])
    updated.priors[in.hpo == T] <- 0.9 / length(official.names[in.hpo == T])
    
    #round updated.priors to avoid floating point error in summing to 1
    priorsequal1 <- 0
    z <- 9
    while(priorsequal1 != 1){
      z <- z +1
      updated.priors <- round(updated.priors, digits = z)
      updated.priors <- updated.priors/sum(updated.priors)
      priorsequal1 <- sum(updated.priors)
    }
    
    hdrda.updated <- hdrda(synd ~ ., data = hdrda.df, prior = updated.priors)
    updated.prediction <- predict(hdrda.updated, newdata = hdrda.df[name.key[names(name.key) %in% input$file1], -1], type = "prob") #what is our pred on the holdout individual given the updated priors 
      
    print(input$hpo)
    print(is.null(input$hpo) | is.null(input$hpo_terms))
    if(is.null(input$hpo) | is.null(input$hpo_terms)) updated.prediction <- predict(hdrda.mod, newdata = hdrda.df[name.key[names(name.key) %in% input$file1], -1], type = "prob")

    return(list(updated.prediction))
  })

  
  output$hpo_synds <- renderDataTable({
    #full hpo list
    print(input$hpo_terms)
    phenotype_2022[phenotype_2022$HPO_ID == names(hpo$name)[which(as.character(hpo$name) == input$hpo_terms)], c(2,5,8)]    #just our dataset HPO
  
    }, options = list(pageLength = 5, autoWidth = F, scrollX = T))
  

  output$diagnosis_scree <- renderPlotly({

    posterior.distribution <- sort(new.mod()[[1]]$posterior, decreasing = T)

    #used to be part of plot.df: ID = as.factor(1:10),
    plot.df <- data.frame(Probs = round(as.numeric(posterior.distribution[1:10]), digits = 4), Syndrome = as.factor(gsub("_", " ", names(posterior.distribution[1:10]))))
    plot.df$Syndrome <- as.character(plot.df$Syndrome)
    plot.df$Syndrome[plot.df$Syndrome == "Control"] <- "Non-syndromic"

    plot_ly(data = plot.df, x = ~Syndrome, y = ~Probs, type = "bar", color = I("grey"), hoverinfo = paste0("Syndrome: ", "x", "<br>", "Probability: ", "y")) %>%
      layout(xaxis = list(tickvals = gsub("_", " ", plot.df$Syndrome), tickangle = 45, ticktext = c(Syndrome = plot.df$Syndrome, Probability = plot.df$Probs), title = "<b>Syndrome</b>"),
             yaxis = list(title = "<b>Class probability</b>"),
             paper_bgcolor='rgba(245, 245, 245, .9)',
             margin = list(b = 125, l = 50, r = 100)
            )


  })

  # #report generator
  # output$report <- downloadHandler(
  #   # For PDF output, change this to "report.pdf"
  #   filename = paste0("FB2_report.html"),
  #   content = function(file) {
  #     # Copy the report file to a temporary directory before processing it, in
  #     # case we don't have write permissions to the current working dir (which
  #     # can happen when deployed).
  #     tempReport <- file.path(tempdir(), "report.Rmd")
  #     file.copy("~/shiny/Classification_demo/report/report.Rmd", tempReport, overwrite = TRUE)
  #     
  #     # Set up parameters to pass to Rmd document
  #     parameters <- list(variables2 = input$scan, lambda = input$lambda, mutant = input$mutant, mag = input$mag)
  #     
  #     # Knit the document, passing in the `params` list, and eval it in a
  #     # child of the global environment (this isolates the code in the document
  #     # from the code in this app).
  #     rmarkdown::render(tempReport, output_file = file,
  #                       params = parameters#,
  #                       # envir = new.env(parent = globalenv())
  #     )
  #   }
  # )
  
  
}


#Run app
shinyApp(ui = ui, server = server)






