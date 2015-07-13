require(ggplot2, quietly = TRUE)
cuts <- seq(0,0.15,by = 0.001)
l <- length(cuts)
n <- 200

mainDir <- getwd()
subDir <- "fdp_output"
k.ind <- 5
#fdp.deseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_deseq2_simseq.RDS"))
fdp.ql <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_simseq_ql.RDS"))
fdp.qlspline <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_simseq_spline.RDS"))
fdp.ql.bart <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_simseq_bart.RDS"))
#fdp.edgeR <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_edger_simseq.RDS"))
#fdp.samseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_samseq_simseq.RDS"))
fdp.voom <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_voom_simseq.RDS"))

fdp1 <- data.frame(fdp = c(rowMeans(fdp.ql),rowMeans(fdp.qlspline),rowMeans(fdp.ql.bart),rowMeans(fdp.voom)),
                   sderr = c(apply(fdp.ql, 1, sd)/sqrt(n),apply(fdp.qlspline, 1, sd)/sqrt(n), apply(fdp.ql.bart, 1, sd)/sqrt(n),
                             apply(fdp.voom, 1, sd)/sqrt(n)),
                   cuts = rep(cuts,4),
                   Method = c(rep("QuasiSeqQL",l),rep("QuasiSeqQLSpline",l),rep("QuasiSeqBart",l), rep("Voom",l)),
                   Simulation = rep("SimSeq", l*4),
                   SS = rep("Sample Size: 5", l*4))

#fdp.deseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_deseq2_nb.RDS"))
fdp.ql <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_nb_ql.RDS"))
fdp.qlspline <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_nb_spline.RDS"))
fdp.ql.bart <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_nb_bart.RDS"))
#fdp.edgeR <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_edger_nb.RDS"))
#fdp.samseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_samseq_nb.RDS"))
fdp.voom <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_voom_nb.RDS"))
fdp2 <- data.frame(fdp = c(rowMeans(fdp.ql),rowMeans(fdp.qlspline),rowMeans(fdp.ql.bart),rowMeans(fdp.voom)),
                   sderr = c(apply(fdp.ql, 1, sd)/sqrt(n),apply(fdp.qlspline, 1, sd)/sqrt(n), apply(fdp.ql.bart, 1, sd)/sqrt(n),
                             apply(fdp.voom, 1, sd)/sqrt(n)),
                   cuts = rep(cuts,4),
                   Method = c(rep("QuasiSeqQL",l),rep("QuasiSeqQLSpline",l),rep("QuasiSeqBart",l), rep("Voom",l)),
                   Simulation = rep("NegBin", l*4),
                   SS = rep("Sample Size: 5", l*4))

k.ind <- 10
#fdp.deseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_deseq2_nb.RDS"))
fdp.ql <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_nb_ql.RDS"))
fdp.qlspline <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_nb_spline.RDS"))
fdp.ql.bart <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_nb_bart.RDS"))
#fdp.edgeR <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_edger_nb.RDS"))
#fdp.samseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_samseq_nb.RDS"))
fdp.voom <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_voom_nb.RDS"))
fdp3 <- data.frame(fdp = c(rowMeans(fdp.ql),rowMeans(fdp.qlspline),rowMeans(fdp.ql.bart),rowMeans(fdp.voom)),
                   sderr = c(apply(fdp.ql, 1, sd)/sqrt(n),apply(fdp.qlspline, 1, sd)/sqrt(n), apply(fdp.ql.bart, 1, sd)/sqrt(n),
                             apply(fdp.voom, 1, sd)/sqrt(n)),
                   cuts = rep(cuts,4),
                   Method = c(rep("QuasiSeqQL",l),rep("QuasiSeqQLSpline",l),rep("QuasiSeqBart",l), rep("Voom",l)),
                   Simulation = rep("NegBin", l*4),
                   SS = rep("Sample Size: 10", l*4))

#fdp.deseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_deseq2_simseq.RDS"))
fdp.ql <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_simseq_ql.RDS"))
fdp.qlspline <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_simseq_spline.RDS"))
fdp.ql.bart <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_simseq_bart.RDS"))
#fdp.edgeR <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_edger_simseq.RDS"))
#fdp.samseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_samseq_simseq.RDS"))
fdp.voom <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_voom_simseq.RDS"))
fdp4 <- data.frame(fdp = c(rowMeans(fdp.ql),rowMeans(fdp.qlspline),rowMeans(fdp.ql.bart),rowMeans(fdp.voom)),
                   sderr = c(apply(fdp.ql, 1, sd)/sqrt(n),apply(fdp.qlspline, 1, sd)/sqrt(n), apply(fdp.ql.bart, 1, sd)/sqrt(n),
                             apply(fdp.voom, 1, sd)/sqrt(n)),
                   cuts = rep(cuts,4),
                   Method = c(rep("QuasiSeqQL",l),rep("QuasiSeqQLSpline",l),rep("QuasiSeqBart",l), rep("Voom",l)),
                   Simulation = rep("SimSeq", l*4),
                   SS = rep("Sample Size: 10", l*4))

k.ind <- 20
#fdp.deseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_deseq2_nb.RDS"))
fdp.ql <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_nb_ql.RDS"))
fdp.qlspline <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_nb_spline.RDS"))
fdp.ql.bart <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_nb_bart.RDS"))
#fdp.edgeR <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_edger_nb.RDS"))
#fdp.samseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_samseq_nb.RDS"))
fdp.voom <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_voom_nb.RDS"))
fdp5 <- data.frame(fdp = c(rowMeans(fdp.ql),rowMeans(fdp.qlspline),rowMeans(fdp.ql.bart),rowMeans(fdp.voom)),
                   sderr = c(apply(fdp.ql, 1, sd)/sqrt(n),apply(fdp.qlspline, 1, sd)/sqrt(n), apply(fdp.ql.bart, 1, sd)/sqrt(n),
                             apply(fdp.voom, 1, sd)/sqrt(n)),
                   cuts = rep(cuts,4),
                   Method = c(rep("QuasiSeqQL",l),rep("QuasiSeqQLSpline",l),rep("QuasiSeqBart",l), rep("Voom",l)),
                   Simulation = rep("NegBin", l*4),
                   SS = rep("Sample Size: 20", l*4))

#fdp.deseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_deseq2_simseq.RDS"))
fdp.ql <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_simseq_ql.RDS"))
fdp.qlspline <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_simseq_spline.RDS"))
fdp.ql.bart <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_simseq_bart.RDS"))
#fdp.edgeR <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_edger_simseq.RDS"))
#fdp.samseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_samseq_simseq.RDS"))
fdp.voom <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_voom_simseq.RDS"))
fdp6 <- data.frame(fdp = c(rowMeans(fdp.ql),rowMeans(fdp.qlspline),rowMeans(fdp.ql.bart),rowMeans(fdp.voom)),
                   sderr = c(apply(fdp.ql, 1, sd)/sqrt(n),apply(fdp.qlspline, 1, sd)/sqrt(n), apply(fdp.ql.bart, 1, sd)/sqrt(n),
                             apply(fdp.voom, 1, sd)/sqrt(n)),
                   cuts = rep(cuts,4),
                   Method = c(rep("QuasiSeqQL",l),rep("QuasiSeqQLSpline",l),rep("QuasiSeqBart",l), rep("Voom",l)),
                   Simulation = rep("SimSeq", l*4),
                   SS = rep("Sample Size: 20", l*4))

fdp <- rbind(fdp1,fdp2,fdp3,fdp4,fdp5,fdp6)
fdp$Method <- as.factor(fdp$Method)
fdp$Simulation <- as.factor(fdp$Simulation)
fdp$SS <- as.factor(fdp$SS)

p <- ggplot(data = fdp) 
p + geom_line(aes(x = cuts, y = fdp, colour = Simulation), size = 0.5, linetype = 11) +
  geom_ribbon(aes(x = cuts, ymax = fdp + 1.96 * sderr, colour = Simulation, fill = Simulation,
                  ymin = fdp - 1.96 * sderr), alpha = 0.3, size = 0.5) +
  scale_y_continuous(breaks = c(-0.05,0,0.05,0.1,0.15)) + 
  scale_x_continuous(breaks = c(0,0.05,0.10)) +
  #   scale_colour_manual(values = c("red", "blue")) +
  #   scale_fill_manual(values = c("red", "blue")) + 
  geom_hline(yintercept = 0, colour = "gold", size = 1, linetype = "dashed") + 
  theme_bw() + facet_grid(SS ~ Method) + 
  theme(axis.text.x = element_text(angle=90)) +
  xlab("Q-value Cutoff") + ylab("(Average FDP) - (Q-Value Cutoff)")
