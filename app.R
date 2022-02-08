#classification demo w/ hpo terms
library(shiny)
library(shinydashboard)
library(Morpho)
library(rgl)
library(ggplot2)
library(caret)
library(sparsediscrim)
library(plotly)
library(shinycssloaders)
library(Rvcg)
library(shinyjs)
library(RvtkStatismo)
library(shinyFiles)
library(mesheR)
library(DT)
options(shiny.maxRequestSize=300*1024^2)


# save(, file = "~/shiny/Classification_demo/demo_objects.Rdata")
# load("/srv/shiny-server/Classification_demo/demo_objects.Rdata")
load("~/app_startup.Rdata")
load("~/adjusted_PCs.Rdata")
# load("~/mshape.Rdata")
atlas$vb[-4,] <- t(synd.mshape)
atlas.lms <- read.mpp("data/atlas_picked_points.pp")
hdrda.mod <- hdrda(synd ~ ., data = hdrda.df[,1:81])

predPC.lm <- function(fit, datamod){
  mat <- model.matrix(datamod)
  pred <- mat %*% fit
  
  return(pred)
}


body <- dashboardBody(
  useShinyjs(),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  fluidRow(
    # box(title = tags$b("About this app"),
    #     solidHeader = T,
    #     status = "warning",
    #     collapsible = T,
    #     collapsed = T,
    #     width = 12,
    #     "This application serves as a demonstration of our syndrome classifier. Choose one of the preloaded faces under the [Pick Parameters] tab. The face will be displayed under the [Scan] tab. The [Syndrome Probabilities] tab will show the top ten syndrome probabilities from the classifier. You can compare the individual's facial shape similarity to the average shape of any syndrome by checking the Make a comparison box under the [Pick parameters] tab."
    # ),
    box(title = tags$b("Scan"),
        width = 6,
        solidHeader = T,
        status = "warning",
        collapsible = F,
        withSpinner(rglwidgetOutput("submitted_face", width = "auto"), type = 6, color = "#757471"),
        hr(style = "border-top: 1.2px solid #d3d3d3;"),
        splitLayout(cellWidths = c("30%", "30%"),
                    numericInput("age", label = "Current age (years)", value = 33, min = .5, max = 80),
                    selectInput("sex", label = "Sex", choices = c("Female", "Male"), selected = "Male")
        ),
        shinyFilesButton("file1", "Upload Files", "Please select a file", multiple = TRUE, viewtype = "detail")
    ),
    box(title = tags$b("Syndrome probabilities"),
        status = "warning",
        solidHeader = T,
        collapsible = F,
        withSpinner(plotlyOutput("diagnosis_scree"), type = 6, color = "#757471")
    ),
    box(title = tags$b("Compare & visualize morphology"),
        width = 6,
        solidHeader = T,
        status = "warning",
        collapsible = F,
        splitLayout(cellWidths = c("50%", "50%"),
                    disabled(checkboxInput("dense", label = "Dense landmarks")),
                    disabled(downloadButton("report", "Download report"))
        ),
        disabled(checkboxGroupInput("compare", label = "", choices = "Make a comparison?")),
        conditionalPanel(condition = "input.compare == 'Make a comparison?'",
                         selectInput("reference", label = "Syndrome", choices = sort(unique(hdrda.df$synd)), selected = "Non-syndromic"), 
                         checkboxInput("displace", label = "Plot displacement vectors?", value = F),
                         conditionalPanel("input.displace == 1",
                                          sliderInput("transparency", label = "Mesh transparency", min = 0, max = 1, step = .1, value = .5)
                         )
        )
    )
  ),
  box(title = tags$b("Human Phenotype Ontology"),
      width = 6,
      status = "warning",
      solidHeader = T,
      collapsible = T,
      collapsed = T,
      checkboxGroupInput("hpo", label = "Extra phenotype info", choices = "Add HPO terms"),
      conditionalPanel(condition = "input.hpo == 'Add HPO terms'",
                       selectInput("hpo_terms", label = "Terms", choices = c("None", as.character(hpo$name[names(hpo$name) %in% unique(phenotype_2022$HPO_ID)])), selected = "None")
      ),
      hr(style = "border-top: 1.2px solid #d3d3d3;"),
      DTOutput('hpo_synds'),
      hr(style = "border-top: 1.2px solid #d3d3d3;"),
      h6("HPO database built from August 2021 release. ",  a("Read more here.", href="https://hpo.jax.org/app/help/introduction"))
  )
)


# We'll save it in a variable `ui` so that we can preview it in the console
dbHeader <- dashboardHeader()
dbHeader$children[[2]]$children <-  tags$a(href='../',
                                           tags$img(src="uc_logo.png", height='30', width='116'))
ui <- dashboardPage(title = "Classification demo",
                    dbHeader,
                    dashboardSidebar(disable = T),
                    body
)

server <- function(input, output, session){
  
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  shinyFileChoose(input, "file1", roots = volumes, session = session)  # by setting `allowDirCreate = FALSE` a user will not be able to create a new directory
  
  mesh.et.lms <- eventReactive(input$file1, {
    parsed.file.path <- parseFilePaths(volumes, input$file1)$datapath
    print(parsed.file.path)
    
    file.mesh <- file2mesh(parsed.file.path[grepl("*.ply", parsed.file.path)])
    # file.name <- substr(input$file1$name, start = 1, stop = nchar(input$file1$name) - 4)
    file.name <- "test"
    
    file.lms <- read.mpp(parsed.file.path[grepl("*.pp", parsed.file.path)])
    
    tmp.fb <- rotmesh.onto(file.mesh, refmat = file.lms, tarmat = atlas.lms, scale = T)
    gp.fb <- tmp.fb$yrot
    
    tmp.fb <- tmp.fb$mesh
    
    Kernels <- IsoKernel(0.1,atlas)
    mymod <- statismoModelFromRepresenter(atlas, kernel = Kernels, ncomp = 100)
    postDef <- posteriorDeform(mymod, tmp.fb, modlm = atlas.lms, samplenum = 1000)
    # return(list(postDef))
    withProgress(message = 'Fitting mesh', value = 0, {
      
      for (i in 1:7) {
        # Increment the progress bar, and update the detail text.
        incProgress(1/7, detail = c(rep("Shaping, thinking, doing...", 4), "Coffee break?", "Non-rigid deformation", "Last few measurements")[i])
        
        if(i < 4) postDef <- posteriorDeform(mymod, tmp.fb, modlm = atlas.lms, tarlm = gp.fb, samplenum = 1000, reference = postDef)
        
        if(i > 4){
          postDefFinal <- postDef
          postDefFinal <- posteriorDeform(mymod, tmp.fb, modlm=atlas.lms, samplenum = 3000, reference = postDefFinal, deform = T, distance = 3)
        }
      }
      
    })
    
    sample1k <- sample(1:27903, 500)
    #register landmarks to the space
    registered.mesh <- rotmesh.onto(postDefFinal, t(postDefFinal$vb[-4, sample1k]), synd.mshape[sample1k,], scale = T)$mesh
    
    return(list(postDefFinal, file.name, file.mesh))
  })
  
  
  output$submitted_face <- renderRglwidget({
    pdf(NULL)
    dev.off()
    
    par3d(userMatrix = diag(4), zoom = .75)
    bg3d(color = rgb(245/255,245/255,245/255,.9))
    plot3d(vcgSmooth(mesh.et.lms()[[1]]), col = "slategrey", axes = F, specular = 1, xlab = "", ylab = "", zlab = "", aspect = "iso")  
    #debug: plot3d(vcgSmooth(jovid), col = "#969FB4", axes = F, specular = 1, xlab = "", ylab = "", zlab = "", aspect = "iso")  
    
    if(input$dense > 0) points3d(t(mesh.et.lms()[[1]]$vb)[,-4], col = 2, alpha = .5)
    
    #     if(length(input$compare) > 0){
    # 
    #       #find individuals' registered lms
    #       ind.lms <- matrix(as.numeric(filtered.lms[as.character(filtered.lms[,1]) %in% mesh.et.lms()[[3]],-1:-2]), nrow = 65, byrow = T)
    # 
    #       #rotonto atlas
    #       scaled.lms <- procOPA(mshape(arrayspecs(filtered.lms[,-1:-2], 65, 3)), as.matrix(atlas.lms))$Bhat
    #       scaled.mesh <- rotmesh.onto(atlas, refmat = as.matrix(atlas.lms), tarmat = as.matrix(scaled.lms), scale = T, reflection = T)
    # 
    #       #calculate specified sydrome mean
    #       syndrome.mean <- matrix(as.numeric(colMeans(filtered.lms[filtered.lms$Syndrome == input$reference,-1:-2])), nrow = 65, byrow = T)
    #       #tps3d
    #       clear3d()
    # 
    #       # par3d(userMatrix = matrix(c(.998,-.005,.0613,0,.0021,.999,.045,0,-.061,-.045,.997,0,0,0,0,1),ncol =4,nrow = 4))
    #       par3d(userMatrix = matrix(c(-.017,-.999,-.022,0,-.999,-.016,-.03,0,.03,-.023,.999,0,0,0,0,1),ncol =4,nrow = 4))
    #       par3d(zoom = .7)
    # 
    # 
    #       synd.mesh <- tps3d(scaled.mesh$mesh, scaled.lms, syndrome.mean)
    #       ind.mesh <- rotmesh.onto(mesh.et.lms()[[1]], as.matrix(mesh.et.lms()[[2]]), ind.lms, scale = T, reflection = T)$mesh
    # 
    #       if(input$displace == F){
    #         bg3d(color = rgb(245/255, 245/255, 245/255, .9))
    #         mD.synd <- meshDist(ind.mesh, synd.mesh , plot = F, scaleramp = F, displace = input$displace, alpha = 1)
    #         a <- render(mD.synd, displace = input$displace, alpha = 1)
    #       } else if(input$displace){
    #         bg3d(color = "#1a1a1a")
    #         mD.synd <- meshDist(ind.mesh, synd.mesh, plot = F, scaleramp = F, displace = input$displace, alpha = input$transparency)
    #         a <- render(mD.synd, displace = input$displace, alpha = input$transparency)
    #       }
    # 
    #     }
    
    rglwidget()
    
  })  
  
  observeEvent(input$file1,
               {
                 enable("dense")
                 enable("compare")
                 enable("report")
               })
  
  #update the priors if given hpo terms
  # new.mod <- eventReactive(input$update_model,{
  new.mod <- reactive({
    #project new face in PC space
    new.scores <- getPCscores(t(mesh.et.lms()[[1]]$vb[-4,]), PC.eigenvectors, synd.mshape)
    # jovid <- file2mesh("data/da_reg.ply")
    
    sample1k <- sample(1:27903, 500)
    #register landmarks to the space
    registered.mesh <- rotmesh.onto(mesh.et.lms()[[1]], t(mesh.et.lms()[[1]]$vb[-4, sample1k]), synd.mshape[sample1k,], scale = T)$mesh
    new.scores <- data.frame(getPCscores(t(registered.mesh$vb[-4,]), PC.eigenvectors, synd.mshape))
    
    #correct submitted face for regression effects
    age.poly <- predict(poly(d.meta$Age,3), newdata = input$age)
    # datamod <- ~ as.numeric("Male" == "Female") + age.poly[1] + age.poly[2] + age.poly[3]
    datamod <- ~ as.numeric(input$sex == "Female") + age.poly[1] + age.poly[2] + age.poly[3]
    expected.values <- predPC.lm(fit = age.sex.lm$coefficients, datamod = datamod)
    
    new.adjusted <- new.scores - expected.values
    colnames(new.adjusted) <- colnames(PC.scores)
    
    if(is.null(input$hpo) | is.null(input$hpo_terms) | input$hpo_terms == "None"){
      updated.prediction <- predict(hdrda.mod, newdata = new.adjusted[1,1:80], type = "prob")
    } else{
      
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
      
      hdrda.updated <- hdrda(synd ~ ., data = hdrda.df[,1:81], prior = updated.priors)
      updated.prediction <- predict(hdrda.updated, newdata = new.adjusted[1:80], type = "prob") #what is our pred on the holdout individual given the updated priors
      
    }
    
    return(list(updated.prediction))
  })
  
  output$hpo_synds <- renderDT({
    #full hpo list
    print(sum(grepl("*.ply", parseFilePaths(volumes, input$file1)$datapath)))
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






