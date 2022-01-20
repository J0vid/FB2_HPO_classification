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

# save(list = ls(), file = "~/shiny/Classification_demo/demo_objects.Rdata")
load("/srv/shiny-server/Classification_demo/demo_objects.Rdata")

#define new morpho plotting method until he releases it on cran

render <- function(x,...) UseMethod("render")

#' @rdname render
#' @method render meshDist
#' @export
render.meshDist <- function(x,from=NULL,to=NULL,steps=NULL,ceiling=NULL,uprange=NULL,tol=NULL,tolcol=NULL,rampcolors=NULL,NAcol=NULL,displace=FALSE,shade=TRUE,sign=NULL,add=FALSE,scaleramp=NULL,...) {
  clost <- x$clost
  dists <- x$dists
  distsOrig <- dists
  colorall <- x$cols
  colramp <- x$colramp
  params <- x$params
  distqual <- x$distqual
  if (!is.null(tolcol))
    tolcol <- colorRampPalette(tolcol)(1)
  if (!add) {
    if (rgl.cur() !=0)
      rgl.clear()
  }
  if (!is.null(from) || !is.null(to) || !is.null(uprange) ||  !is.null(tol)  ||  !is.null(sign) || !is.null(steps) || !is.null(rampcolors) || !is.null(NAcol) || !is.null(tolcol) || !is.null(scaleramp)) {
    neg=FALSE
    colMesh <- x$colMesh
    if(is.null(steps))
      steps <- x$params$steps
    if (is.null(rampcolors))
      rampcolors <- x$params$rampcolors
    if (is.null(NAcol))
      NAcol <- x$params$NAcol
    if (is.null(tolcol))
      tolcol <- x$params$tolcol
    if (is.null(tol))
      tol <- x$params$tol
    if(is.null(sign))
      sign <- x$params$sign
    if (!sign) {
      distsOrig <- dists
      dists <- abs(dists)
    }
    if(is.null(ceiling))
      ceiling <- x$params$ceiling
    if(is.null(uprange))
      uprange <- x$params$uprange
    
    if (is.null(from)) {
      mindist <- min(dists)
      if (sign && mindist < 0 ) {
        from <- quantile(dists,probs=(1-uprange)) 
        neg <- TRUE            
      } else {
        from <- 0
      }             
    }
    if (is.null(scaleramp))
      scaleramp <- x$params$scaleramp
    
    if (from < 0)
      neg <- TRUE
    if (is.null(to))
      to <- quantile(dists,probs=uprange)    
    if(ceiling)
      to <- ceiling(to)
    
    to <- to+1e-10
    #ramp <- blue2green2red(maxseq*2)
    ramp <- colorRampPalette(rampcolors)(steps-1)
    colseq <- seq(from=from,to=to,length.out=steps)
    coldif <- colseq[2]-colseq[1]
    if (neg && sign) {
      
      negseq <- length(which(colseq<0))
      poseq <- steps-negseq
      maxseq <- max(c(negseq,poseq))
      if (scaleramp) {
        ramp <- colorRampPalette(rampcolors)(maxseq*2)
        ramp <- ramp[c(maxseq-negseq+1):(maxseq+poseq)]
        
      }
      else
        ramp <- colorRampPalette(rampcolors)(steps-1)
      distqual <- ceiling(((dists+abs(from))/coldif)+1e-14)
      #distqual[which(distqual < 1)] <- steps+10
    } else if (from > 0) {
      distqual <- ceiling(((dists-from)/coldif)+1e-14)
    } else {
      distqual <- ceiling((dists/coldif)+1e-14)
    }
    distqual[which(distqual < 1)] <- steps+10
    colorall <- ramp[distqual]
    if (!is.null(tol)) {
      if ( length(tol) < 2 ) {
        if (sign) {
          tol <- c(-tol,tol)
        } else {
          tol <- c(0,tol)
        }
      }
      good <- which(abs(dists) < tol[2])
      colorall[good] <- tolcol
    }
    colfun <- function(x){x <- colorall[x];return(x)}
    colMesh$material$color <- matrix(colfun(colMesh$it),dim(colMesh$it))
    colMesh$material$color[is.na(colMesh$material$color)] <- NAcol
    #colMesh$material$color <- matrix(colfun(colMesh$it),dim(colMesh$it))
    colramp <- list(1,colseq, matrix(data=colseq, ncol=length(colseq),nrow=1),col=ramp,useRaster=T,ylab="Distance in mm",xlab="",xaxt="n")
  } else {
    if (is.null(tol))
      tol <- x$params$tol
    colramp <- x$colramp
    colMesh <- x$colMesh
  }
  if (is.null(tolcol))
    tolcol <- x$params$tolcol
  
  if (shade)
    shade3d(vcgUpdateNormals(colMesh),specular="black",meshColor="legacy",...)
  if (displace) {
    dismesh <- colMesh
    vl <- dim(colMesh$vb)[2]
    dismesh$vb <- cbind(colMesh$vb,rbind(clost,1))
    dismesh$it <- rbind(1:vl,1:vl,(1:vl)+vl)
    dismesh$material$color <- rbind(colorall,colorall,colorall)
    wire3d(dismesh,lit=FALSE,meshColor="legacy")
  }
  diffo <- ((colramp[[2]][2]-colramp[[2]][1])/2)
  image(colramp[[1]],colramp[[2]][-1]-diffo,t(colramp[[3]][1,-1])-diffo,col=colramp[[4]],useRaster=TRUE,ylab="Distance in mm",xlab="",xaxt="n")
  if (!is.null(tol)) {
    if (sum(abs(tol)) != 0)
      image(colramp[[1]],c(tol[1],tol[2]),matrix(c(tol[1],tol[2]),1,1),col=tolcol,useRaster=TRUE,add=TRUE)
  }
  params <- list(steps=steps,from=from,to=to,uprange=uprange,ceiling=ceiling,sign=sign,tol=tol,rampcolors=rampcolors,NAcol=NAcol,tolcol=tolcol)
  out <- list(colMesh=colMesh,dists=distsOrig,cols=colorall,colramp=colramp,params=params,distqual=distqual,clost=clost)
  
  class(out) <- "meshDist"
  invisible(out)
}


body <- dashboardBody(
  tags$style(HTML("

body {
background: rgba(245, 245, 245, .9);
}

.skin-blue .main-header .logo:hover {
background-color: #901628;
}

.box.box-solid.box-warning>.box-header {
  color:#fff;
  background:#a6192e
                    }

.box.box-solid.box-warning{
background: rgba(245, 245, 245, .9);
border-bottom-color:#666666;
border-left-color:#666666;
border-right-color:#666666;
border-top-color:#666666;
box-shadow: 0 4px 8px 0 rgba(0, 0, 0, 0.9);
border-top-left-radius:5px;
border-top-right-radius:5px;
}

.box.box-solid.box-warning:hover {
    box-shadow: 0 12px 24px 0 rgba(0, 0, 0, 0.9);
    transition: box-shadow 0.3s ease-in-out;
}

.skin-blue .main-header .navbar {
    background-color: #a6192e;
}

.skin-blue .main-header .logo {
    background-color: #a6192e;
    color: #fff;
    border-bottom: 0 solid transparent;
}

                                    ")),
  fluidRow(
    box(title = tags$b("About this app"),
        solidHeader = T,
        status = "warning",
        collapsible = T,
        width = 12,
        #verbatimTextOutput("about"),
        "This application serves as a demonstration of our syndrome classifier. Choose one of the preloaded faces under the [Pick Parameters] tab. The face will be displayed under the [Scan] tab along with the 65 3D landmarks. The [Syndrome Probabilities] tab will show the top ten syndrome probabilities from the classifier. You can compare the individual's facial shape similarity to the average shape of any syndrome by checking the Make a comparison box under the [Pick parameters] tab."
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
        selectInput("file1", "", choices = c("Apert", "Achondroplasia", "Williams", "Treacher Collins", "Noonan", "22q deletion", "Cockayne", "Nager"), selected = "22q deletion", multiple = FALSE, selectize = TRUE, width = NULL, size = NULL),
        splitLayout(cellWidths = c("50%", "50%"),
                    checkboxGroupInput("compare", label = "Morphology", choices = "Make a comparison?"),
                    checkboxInput("dense", label = "Dense landmarks")
                    ),
        splitLayout(cellWidths = c("50%", "50%"),
                    checkboxGroupInput("hpo", label = "Extra phenotype info", choices = "Add HPO terms"),
                    disabled(downloadButton("report", "Send me the results!"))
                    ),
        conditionalPanel(condition = "input.compare == 'Make a comparison?'",
                         selectInput("reference", label = "Syndrome", choices = sort(unique(FB2_class$Syndrome)), selected = "Non-syndromic"), 
                         checkboxInput("displace", label = "Plot displacement vectors?", value = F),
                         conditionalPanel("input.displace == 1",
                                          sliderInput("transparency", label = "Mesh transparency", min = 0, max = 1, step = .1, value = .5)
                         )
        ),
        conditionalPanel(condition = "input.hpo == 'Add HPO terms'",
                         textInput("hpo_input", label = "HPO", value = "limb"),
                         uiOutput('variables')
                         # actionButton("update_model", "Update model")
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
  
  # output$about <- renderText({"The following app is a demonstration of our syndrome classifier with a selection of faces representing different genetic syndromes. When a face is selected, it will be shown in the scan window along with the 65 3D points we measure. The facial measurements will be run through the syndrome classifier and the top 10 syndromic probabilities will be presented."})
  
  name.key <- c(Apert = "140903101436.ply", Achondroplasia = "150706091510.ply", Williams = "150715095947.ply", `Treacher Collins` = "150910091215.ply", Noonan = "160722120720.ply", `22q deletion` = "160728105214.ply", Cockayne = "170728172355.ply", Nager = "170823151810.ply")

  outVar <- reactive({
    hpo.terms <- hpo$name[grep(tolower(hpo$name), pattern = tolower(input$hpo_input))]
    hpo.id <- names(hpo.terms)
    hpo.terms <- as.character(hpo.terms)
    return(list(hpo.id, hpo.terms))
  })
  
  #slow down reactivity of text input
  ontology.list <- debounce(outVar, 1200)
  
  output$variables <- renderUI({
    selectInput('variables2', 'Associated HPO terms', ontology.list()[[2]], multiple = T)
  })
  
  
  mesh.et.lms <- eventReactive(input$file1, {
  # mesh.et.lms <- reactive({
    # file1 <- name.key[names(name.key) %in% "Nager"]
    file1 <- name.key[names(name.key) %in% input$file1]
    file.mesh <- file2mesh(paste0("/srv/shiny-server/Classification_demo/scans/", file1))

    file.name <- substr(file1, start = 1, stop = nchar(file1) - 4)
    file.lms <- read.table(paste0("/srv/shiny-server/Classification_demo/lms/", file.name, ".txt"))
    
    return(list(file.mesh, file.lms, file.name))
  })

  output$LM_face <- renderRglwidget({
    #look at selected person's file name, find the scan and paste it here
    scan.name <- substr(FB2_class$Subj_ID2, 1, stop = nchar(as.character(FB2_class$Subj_ID2))-2)
    
    clear3d()
    bg3d(color = rgb(245/255, 245/255, 245/255, .9))
    r.at.lms <- procOPA(as.matrix(atlas.lms), as.matrix(mesh.et.lms()[[2]]))$R
    r.at.lms <- as.matrix(mesh.et.lms()[[2]]) %*% as.matrix(r.at.lms)
    # par3d(userMatrix = matrix(c(.998,-.005,.0613,0,.0021,.999,.045,0,-.061,-.045,.997,0,0,0,0,1),ncol =4,nrow = 4))
    par3d(userMatrix = diag(4))
    aspect3d("iso")
    par3d(zoom = .7)
    
    meshrot <- rotmesh.onto(mesh.et.lms()[[1]], as.matrix(mesh.et.lms()[[2]]), as.matrix(atlas.lms))
    
    plot.scan <- plot3d(meshrot$mesh, col = adjustcolor("lightgrey", .3), alpha = .3, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "")
    spheres3d(as.matrix(meshrot$yrot), radius = .75, color = "red")
    
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
  
  #change to enable when a scan is uploaded: DONE
  observeEvent(is.null(input$file1) == F,
               {
                 disable("report")
               })
  
  #update the priors if given hpo terms
  # new.mod <- eventReactive(input$update_model,{
  new.mod <- reactive({
    updated.priors <- table(filtered.lms$Syndrome)
    
    #zero out known impossibilities based on HPO
    #when someone selects specific terms, find the omims/syndrome name in hpo_pos
    #here's an example with facial asymmetry selected: hpo_pos[hpo_pos[,5] == "HP:0000322",]
    print(hpo_pos[hpo_pos[,5] == names(hpo$name)[hpo$name == input$variables2],3])
    
    in.hpo <- rep(NA, length(unique(filtered.lms$Syndrome)))
    
    
    
    # View(data.frame(updated.priors, in.hpo[sort(official.names, index.return = T)$ix]))
    
    #resort to match the table and select the right subset
    # in.hpo <- in.hpo[sort(official.names, index.return = T)$ix]
    # 
    # updated.priors[in.hpo == F] <- 0
    # #scale priors
    # updated.priors <- updated.priors/sum(updated.priors)
    
    cva.data <- data.frame(ID = filtered.lms$ID, Syndrome = as.character(filtered.lms$Syndrome), scores = bad.cva$CVscores)
    mod.data <- data.frame(Syndrome = as.character(filtered.lms$Syndrome), scores = bad.cva$CVscores)
    
    official.names <- levels(mod.data$Syndrome)
    official.names[10] <- "CLEFT LIP/PALATE"
    filtered.lms$Syndrome[filtered.lms$Syndrome == "Cleft_Lip_Palate"] <- "CLEFT LIP/PALATE"
    official.names[2] <- "CHROMOSOME 1P36 DELETION SYNDROME"
    official.names[3] <- "22q11.2 DELETION SYNDROME"
    official.names[4] <- "CRI-DU-CHAT SYNDROME"
    official.names[5] <- "ACHONDROPLASIA"
    official.names[6] <- "ANGELMAN SYNDROME"
    official.names[7] <- "APERT SYNDROME"
    official.names[8] <- "CARDIOFACIOCUTANEOUS SYNDROME"
    official.names[9] <- "CHARGE SYNDROME"
    official.names[10] <- "CLEFT LIP/PALATE"
    official.names[11] <- "COCKAYNE SYNDROME"
    official.names[12] <- "COFFIN-SIRIS SYNDROME"
    official.names[13] <- "COHEN SYNDROME"
    official.names[14] <- "CORNELIA DE LANGE SYNDROME"
    official.names[15] <- "COSTELLO SYNDROME"
    official.names[16] <- "CROUZON SYNDROME"
    official.names[17] <- "DOWN SYNDROMETRISOMY 21"
    official.names[18] <- "EHLERS-DANLOS SYNDROME"
    official.names[19] <- "FRAGILE X"
    official.names[21] <- "JACOBSEN SYNDROME"
    official.names[22] <- "JOUBERT SYNDROME"
    official.names[23] <- "KABUKI SYNDROME"
    official.names[25] <- "LOEYS-DIETZ SYNDROME"
    official.names[26] <- "MARFAN SYNDROME; MFS"
    official.names[27] <- "MOEBIUS SYNDROME; MBS"
    official.names[28] <- "MUCOPOLYSACCHARIDOSIS"
    official.names[29] <- "ACROFACIAL DYSOSTOSIS 1, NAGER TYPE; AFD1"
    official.names[30] <- "NEUROFIBROMATOSIS"
    official.names[32] <- "NOONAN SYNDROME"
    official.names[33] <- "OSTEOGENESIS IMPERFECTA"
    official.names[34] <- "PHELAN-MCDERMID SYNDROME"
    official.names[35] <- "PIERRE ROBIN"
    official.names[36] <- "PITT-HOPKINS"
    official.names[37] <- "PSEUDOACHONDROPLASIA"
    official.names[38] <- "RETT SYNDROME"
    official.names[39] <- "RHIZOMELIC CHONDRODYSPLASIA PUNCTATA"
    official.names[40] <- "RUBINSTEIN-TAYBI SYNDROME"
    official.names[41] <- "SILVER-RUSSELL SYNDROME"
    official.names[42] <- "SMITH-LEMLI-OPITZ SYNDROME"
    official.names[43] <- "SMITH-MAGENIS SYNDROME"
    official.names[44] <- "SOTOS SYNDROME"
    official.names[45] <- "SPONDYLOEPIPHYSEAL"
    official.names[46] <- "STICKLER SYNDROME"
    official.names[47] <- "TREACHER COLLINS"
    official.names[48] <- "TRISOMY 18"
    official.names[49] <- "TURNER"
    official.names[50] <- "VAN DER WOUDE SYNDROME"
    official.names[51] <- "WILLIAMS-BEUREN SYNDROME"
    official.names[52] <- "ECTODERMAL DYSPLASIA 1, HYPOHIDROTIC, X-LINKED"
    
    
    
    for(i in 1:length(in.hpo)) in.hpo[i] <- length(grep(official.names[i], x = hpo_pos[hpo_pos[,5] == names(hpo$name)[hpo$name == input$variables2],3])) > 0
    
    in.hpo[official.names == "Non-syndromic"] <- TRUE
    
    # #old code that zeroed out priors
    # mod.data <- cva.data[filtered.lms$Syndrome %in% official.names[in.hpo == T],-1]
    # mod.data$Syndrome <- factor(mod.data$Syndrome, levels = official.names[in.hpo == T])
    
    #priors are adjusted to uniformly sharing 20% for each syndrome not in the selected HPO, the rest of the weight is uniform with the remaining syndromes in the HPO list

    updated.priors <- rep(NA, length(official.names))
    names(updated.priors) <- official.names
    
    updated.priors[in.hpo == F] <- .05 / length(official.names[in.hpo == F])
    updated.priors[in.hpo == T] <- .95 / length(official.names[in.hpo == T])
    
    updated.priors <- as.numeric(scale(updated.priors, center = F, scale = 1))
    
    hpo.dt <- hpo_pos[hpo_pos[,5] == names(hpo$name)[hpo$name == input$variables2], c(3,5)]
   #currently bugged because it's using exact matches
   #  hpo.dt[,2] <- hpo_pos[hpo_pos[,5] == names(hpo$name)[hpo$name == input$variables2],3] %in% official.names[in.hpo == T]
   # colnames(hpo.dt) <- c("Syndrome name", "In classification data?", "HPO term")
    colnames(hpo.dt) <- c("Syndrome name", "HPO term")
    
    
   
    print(updated.priors)
    
    hdrda.updated <- hdrda(Syndrome ~ ., data = mod.data, prior = updated.priors)
    
    updated.prediction <- predict(hdrda.updated, newdata = rbind(cva.data[cva.data$ID == mesh.et.lms()[[3]],-1:-2], cva.data[cva.data$ID == mesh.et.lms()[[3]],-1:-2]), type = "prob")$post[1,]
    
    if(is.null(input$hpo) | is.null(input$variables2)) updated.prediction <- predict(hdrda.cva, newdata = rbind(cva.data[cva.data$ID == mesh.et.lms()[[3]],-1:-2], cva.data[cva.data$ID == mesh.et.lms()[[3]], -1:-2]), type = "prob")[1,] 
    
    return(list(updated.prediction, hpo.dt, official.names[in.hpo == T]))
  })
  
  output$hpo_synds <- renderDataTable({
    #full hpo list
    new.mod()[[2]]
    #just our dataset HPO
    # new.mod()[[3]]
    
    }, options = list(pageLength = 5))
  
  output$diagnosis_scree <- renderPlotly({
    
    cva.data <- data.frame(ID = filtered.lms$ID, Syndrome = filtered.lms$Syndrome, scores = bad.cva$CVscores)
    #classify individual's scores using the model
    
    # if(input$hpo == F) posterior.distribution <- predict(hdrda.cva, newdata = rbind(cva.data[cva.data$ID == mesh.et.lms()[[3]],-1:-2],cva.data[cva.data$ID == mesh.et.lms()[[3]],-1:-2]), type = "prob")[1,]
    #for testing: posterior.distribution <- predict(hdrda.cva, newdata = rbind(cva.data[cva.data$ID == "150715095947",-1:-2],cva.data[cva.data$ID == "150715095947",-1:-2]), type = "prob")[1,]
    
    posterior.distribution <- new.mod()[[1]]
    
    posterior.distribution <- sort(posterior.distribution, decreasing = T)
    
    #used to be part of plot.df: ID = as.factor(1:10), 
    plot.df <- data.frame(Probs = round(as.numeric(posterior.distribution[1:10]), digits = 4), Syndrome = as.factor(gsub("_", " ", names(posterior.distribution[1:10]))))
    plot.df$Syndrome <- as.character(plot.df$Syndrome)
    plot.df$Syndrome[plot.df$Syndrome == "Control"] <- "Non-syndromic"
    
      # p <- ggplot(data = plot.df) +
      # geom_bar(stat = "identity", aes(x = Syndrome, y = Probs)) +
      # theme_bw() +
      # xlab("Syndrome") +
      # ylab("Class probability")  
      # # scale_x_discrete(labels = gsub("_", " ", plot.df$Syndrome)) 
      # 
      # 
      # 
      # ggplotly(p + 
      #           theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 15),
      #                    axis.text.y = element_text(size = 14),
      #                    axis.title = element_text(size = 15, face = "bold"),
      #                    plot.background = element_rect(fill = rgb(245/255,245/255,245/255,.9), colour = rgb(245/255,245/255,245/255,.9))),
      #          tooltip = c("Syndrome", "Probs"))
      # 
     
      plot_ly(data = plot.df, x = ~Syndrome, y = ~Probs, type = "bar", color = I("grey"), hoverinfo = paste0("Syndrome: ", "x", "<br>", "Probability: ", "y")) %>%
        layout(xaxis = list(tickvals = gsub("_", " ", plot.df$Syndrome), tickangle = 45, ticktext = c(Syndrome = plot.df$Syndrome, Probability = plot.df$Probs), title = "<b>Syndrome</b>"),
               yaxis = list(title = "<b>Class probability</b>"),
               paper_bgcolor='rgba(245, 245, 245, .9)',
               margin = list(b = 125, l = 50, r = 100)
               )
      
    
  })
  
  #report generator
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = paste0("FB2_report.html"),
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("~/shiny/Classification_demo/report/report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      parameters <- list(variables2 = input$scan, lambda = input$lambda, mutant = input$mutant, mag = input$mag)
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = parameters#,
                        # envir = new.env(parent = globalenv())
      )
    }
  )
  
  
}


#Run app
shinyApp(ui = ui, server = server)






