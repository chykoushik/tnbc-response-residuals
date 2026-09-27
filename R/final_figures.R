figure09 <- function(number, draw, width=12, height=9) {
  for(ext in c("pdf","png")) {
    f<-path("figures/final",paste0("Figure",number,".",ext))
    if(ext=="pdf") pdf(f,width=width,height=height,useDingbats=FALSE) else
      png(f,width=width,height=height,units="in",res=300,type="cairo")
    tryCatch({par(family="sans",cex=.85,mar=c(4,4,3,1));draw()},finally=dev.off())
  }
}
figures09 <- function() {
  labels<-c("Immune / IFN","Proliferation","Stromal / EMT","DNA repair","Hypoxia / stress","Apoptosis")
  cols<-c("#0072B2","#009E73","#D55E00","#CC79A7","#E69F00","#56B4E9")
  figure09(1,function() {
    par(mar=c(1,1,3,1));plot.new();plot.window(xlim=c(0,10),ylim=c(0,10))
    title("Study design and frozen analysis scope",font.main=2)
    box<-function(x,y,w,h,textvalue,color="#EAF2F8") {
      rect(x-w/2,y-h/2,x+w/2,y+h/2,col=color,border="#627080")
      text(x,y,textvalue,cex=.95)
    }
    box(2.5,8.8,4.3,1.3,"Public NeoTRIP GSE319641\n401 samples; 251 patients\n241 baseline + 160 D1C2")
    box(7.5,8.8,4.3,1.3,"Public-data and identity audits\nComposition reproduced with InstaPrism\nSix frozen Hallmark programmes")
    box(2.5,6.4,4.3,1.6,"150 complete baseline-D1C2 pairs\n79 CT / 71 CT + atezolizumab\n63 pCR / 87 residual disease")
    box(7.5,6.4,4.3,1.6,"Expected change: nested ridge regression\nArm + baseline pathways + composition\n10 outer folds; outcomes excluded")
    box(5,4.1,8.4,1.3,"Residual = observed change - cross-fitted expected change\nEach patient excluded from its expected-response model")
    box(2.5,1.9,4.3,1.6,"Primary: six adjusted residual-pCR models\nBonferroni family; zero confirmed pathways\nRobustness: hypothesis-generating only")
    box(7.5,1.9,4.3,1.6,"Secondary: repeated nested prediction\nOrthogonal single-cell interpretation\nExternal transportability: SKIPPED")
    arrows(2.5,8.1,2.5,7.3,length=.1);arrows(7.5,8.1,7.5,7.3,length=.1)
    arrows(2.5,5.5,4,4.9,length=.1);arrows(7.5,5.5,6,4.9,length=.1)
    arrows(4,3.4,2.5,2.8,length=.1);arrows(6,3.4,7.5,2.8,length=.1)
    text(5,.35,"Single-cohort inference; no clinical utility or external-validation claim",cex=.9)
  })
  pair<-read09("data/processed/pathways/neotrip_confirmatory_pair_scores_v1.tsv")
  comp<-read09("data/processed/composition/neotrip_pair_composition_change.tsv")
  figure09(2,function() {
    par(mfrow=c(2,2),mar=c(8,5,3,1),mgp=c(2.5,.8,0))
    for(prefix in c("baseline_","d1c2_","delta_")) {
      values<-lapply(axes07,function(a)pair[[paste0(prefix,a)]])
      boxplot(values,names=labels,las=2,col=adjustcolor(cols,.5),outline=FALSE,
        ylab="Frozen singscore",main=switch(prefix,baseline_="A  Baseline",d1c2_="B  D1C2",delta_="C  Observed change"))
      if(prefix=="delta_")abline(h=0,lty=2,col="grey50")
    }
    values<-lapply(c("cancer_epithelial","immune_fraction","stromal_fraction"),
      function(c)comp[[paste0("delta_",c)]])
    boxplot(values,names=c("Cancer epithelial","Immune","Stromal"),las=2,col="#C7DDE8",
      outline=FALSE,ylab="Change (percentage points)",main="D  Composition changes")
    abline(h=0,lty=2,col="grey50")
  })
  residual<-read09("data/processed/residuals/neotrip_crossfitted_residuals_long_v1.tsv")
  figure09(3,function() {
    par(mfrow=c(3,4),mar=c(4.5,4,2,1),mgp=c(2.3,.7,0),cex=.7)
    for(a in axes07) {
      z<-residual[residual$axis==a,];i<-match(a,axes07)
      lim<-range(c(z$observed_delta,z$expected_delta))
      plot(z$expected_delta,z$observed_delta,pch=16,col=adjustcolor(cols[i],.5),
        xlim=lim,ylim=lim,xlab="Expected change (OOF)",ylab="Observed change",main=labels[i])
      abline(0,1,lty=2,col="grey40")
      hist(z$residual,breaks=20,col=adjustcolor(cols[i],.6),border="white",
        xlab="Residual change",main="Residual distribution");abline(v=0,lty=2)
    }
  },width=14,height=10)
  primary<-read09("results/stage05_primary/primary_residual_associations_v2.tsv")
  adjustment<-read09("results/stage05_primary/composition_adjustment_comparison_v2.tsv")
  forest<-function(z,main,color) {
    z<-z[match(axes07,z$axis),]; y<-6:1
    plot(z$odds_ratio,y,log="x",xlim=range(c(z$ci95_lower,z$ci95_upper,.5,2)),ylim=c(.5,6.5),
      yaxt="n",xlab="Odds ratio per residual SD (95% CI)",ylab="",pch=16,col=color,main=main)
    axis(2,y,labels,las=1,cex.axis=.8);segments(z$ci95_lower,y,z$ci95_upper,y,col=color,lwd=2);abline(v=1,lty=2)
  }
  figure09(4,function() {
    par(mfrow=c(1,2),mar=c(5,10,4,2))
    forest(primary,"A  Composition-adjusted primary models","#0072B2")
    mtext("Zero of six pass Bonferroni; 95% intervals are unadjusted",side=3,line=.4,cex=.65)
    base<-adjustment[adjustment$model=="base_model",]
    forest(base,"B  Models before composition adjustment","#D55E00")
    mtext("For comparison only; primary inference remains adjusted",side=3,line=.4,cex=.65)
  },width=15,height=6)
  sensitivity<-read09("results/stage06_robustness/five_primary_sensitivities.tsv")
  permutation<-read09("results/stage06_robustness/permutation_summary.tsv")
  influence<-read09("results/stage06_robustness/influence_summary.tsv")
  random<-read09("results/stage06_robustness/random_gene_null_summary.tsv")
  figure09(5,function() {
    par(mfrow=c(2,2),mar=c(5,9,3,1))
    plot(c(min(sensitivity$beta),max(sensitivity$beta)+.2),c(.5,6.5),type="n",yaxt="n",xlab="Residual association coefficient",ylab="",
      main="A  Five prespecified sensitivities");axis(2,6:1,labels,las=1,cex.axis=.75);abline(v=0,lty=2)
    ids<-unique(sensitivity$sensitivity_id)
    for(i in seq_along(ids)) {
      z<-sensitivity[sensitivity$sensitivity_id==ids[i],];z<-z[match(axes07,z$axis),]
      points(z$beta,6:1+(i-3)*.1,pch=14+i,col=cols[i])
    }
    legend("bottomright",sub("_.*","",ids),col=cols[1:5],pch=15:19,bty="n",cex=.65)
    for(kind in c("permutation","influence","random")) {
      z<-switch(kind,permutation=permutation,influence=influence,random=random);z<-z[match(axes07,z$axis),]
      v<-switch(kind,permutation=z$permutation_p,influence=z$max_relative_influence,random=z$primary_percentile)
      plot(v,6:1,yaxt="n",ylim=c(.5,6.5),pch=16,col=cols,xlab=switch(kind,
        permutation="Frozen permutation p (10,000 per pathway)",influence="Maximum relative coefficient influence",
        random="Matched random-gene percentile"),ylab="",main=switch(kind,permutation="B  Permutation results",
        influence="C  Leave-one-patient-out influence",random="D  Matched random-gene controls"))
      axis(2,6:1,labels,las=1,cex.axis=.75)
    }
  },width=13,height=10)
  if (!file.exists(path("results/stage09_final/single_cell_completion.tsv"))) return(invisible(NULL))
  prediction<-read09("results/stage07_prediction/repeat_metrics.tsv")
  complete<-read09("results/stage09_final/single_cell_completion.tsv")$status=="COMPLETED"
  figure09(6,function() {
    par(mfrow=c(2,2),mar=c(5,5,3,1))
    boxplot(roc_auc~model,prediction,col="#C7DDE8",outline=FALSE,ylim=range(prediction$roc_auc),xlab="Model",ylab="OOF ROC-AUC",
      main="A  Secondary prediction: 10 repeats")
    for(i in 1:5)points(rep(i,10),prediction$roc_auc[prediction$model==paste0("P",i-1)],pch=16,col="#0072B2",cex=.65)
    if(complete) {
      interpretation<-read09("results/stage09_final/pathway_compartment_interpretation.tsv")
      for(time in c("baseline","post")) {
        matrix<-sapply(c("immune","stromal","epithelial","other"),function(comp) {
          z<-interpretation[interpretation$broad_compartment==comp,];z[[paste0(time,"_median_across_patients")]][match(axes07,z$axis)]
        })
        par(mar=c(5,9,3,1))
        image(1:4,1:6,t(matrix[6:1,]),axes=FALSE,col=hcl.colors(30,"YlOrRd",rev=TRUE),
          zlim=c(0,1),xlab="",ylab="",main=if(time=="baseline")"B  Pre: median UCell" else "C  Post: median UCell")
        axis(1,1:4,c("Immune","Stromal","Epithelial","Other"),cex.axis=.7)
        axis(2,1:6,rev(labels),las=1,cex.axis=.7)
        mtext("Three patients; patient-compartment summaries require >=50 cells",side=3,line=.3,cex=.6)
        for(i in 1:4)for(j in 1:6)if(is.finite(matrix[7-j,i]))text(i,j,sprintf("%.2f",matrix[7-j,i]),cex=.75)
      }
    } else {
      plot.new();text(.5,.5,"B  Single-cell support omitted\nTechnical reason recorded",cex=1)
      plot.new();text(.5,.5,"C  No external-validation panel\nGSE122630: SKIPPED",cex=1)
    }
    par(mar=c(2,2,3,1));plot.new();title("D  Final evidence boundaries")
    text(0, .9, "150 NeoTRIP patients; six prespecified pathways\n\n0 multiplicity-confirmed pathways\n\n4 hypothesis-generating; 2 not robust\n\nNo consistent incremental predictive gain\n\nSingle-cell: descriptive only; scoring discordance retained\n\nNo clinical utility or independent replication claim",
      adj=c(0,1),cex=.85)
  },width=13,height=10)
  write09(data.frame(figure=paste0("Figure",1:6),scope=c(
    "Study design, public audit, patient flow and residual construction","Baseline, D1C2, pathway and composition changes",
    "Cross-fitted expected versus observed changes and residual distributions","Primary residual associations and composition adjustment comparison",
    "Five sensitivities, permutations, influence and matched random controls",
    if(complete)"Secondary prediction and orthogonal single-cell interpretation" else "Secondary prediction and evidence synthesis"),
    raster_dpi=300,vector_pdf=TRUE,significance_stars=FALSE),
    "metadata/final_gate/final_figure_plan_v1.tsv")
}
