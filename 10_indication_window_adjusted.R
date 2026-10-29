# Standalone indication/common-window active comparators and Firth adjusted ROR.
# Run source('10_indication_window_adjusted.R'); run_indication_window_adjusted()
# See ANALYSIS_PLAN.txt and VALIDATION.md in the final output.
load_one_rdata <- function(path, table_name) {
  e <- new.env(parent = baseenv()); loaded <- load(path, envir = e)
  choices <- loaded[vapply(loaded, function(n) is.data.frame(e[[n]]), logical(1))]
  selected <- if (table_name %in% choices) table_name else choices
  stopifnot(length(selected) == 1L)
  table <- data.table::as.data.table(e[[selected]])
  data.table::setnames(table, toupper(names(table))); table
}
indication_source_quarter <- function(value) {
  token <- regmatches(toupper(as.character(value)), regexec('([0-9]{2,4})Q([1-4])', toupper(as.character(value))))
  vapply(token, function(x) {
    if (length(x) != 3L) return(NA_integer_)
    year <- as.integer(x[2]); if (year < 100) year <- year + 2000L
    if (year < 2004 || year > 2099) return(NA_integer_)
    year * 4L + as.integer(x[3]) - 1L
  }, integer(1))
}
indication_class_for_key <- function(exposure_keys, indications, diabetes_terms, weight_terms) {
  keys <- unique(data.table::copy(exposure_keys))
  keys[, DRUG_SEQ := trimws(as.character(DRUG_SEQ))]
  indi <- data.table::copy(indications)
  indi[, DRUG_SEQ := trimws(as.character(DRUG_SEQ))]
  unknown_terms <- c('', 'UNKNOWN', 'DRUG USE FOR UNKNOWN INDICATION', 'PRODUCT USED FOR UNKNOWN INDICATION', 'UNKNOWN INDICATION', 'NOT REPORTED', 'NOT SPECIFIED')
  indi[, `:=`(diabetes = INDI_PT %chin% diabetes_terms, weight = INDI_PT %chin% weight_terms,
    specified = !is.na(INDI_PT) & !INDI_PT %chin% unknown_terms)]
  flags <- indi[!is.na(DRUG_SEQ) & nzchar(DRUG_SEQ), .(diabetes = any(diabetes), weight = any(weight),
    specified = any(specified), linked_indi_rows = .N), by = .(PRIMARYID, DRUG_SEQ)]
  keys <- merge(keys, flags, by = c('PRIMARYID', 'DRUG_SEQ'), all.x = TRUE)
  keys[is.na(linked_indi_rows), `:=`(diabetes = FALSE, weight = FALSE, specified = FALSE, linked_indi_rows = 0L)]
  # Unknown means no specified indication, NOT a known other indication.
  keys[, indication_status := data.table::fcase(diabetes & weight, 'both', diabetes, 'diabetes', weight, 'weight_management',
    specified, 'other_known', default = 'unknown')]
  data.table::rbindlist(list(keys[diabetes == TRUE][, indication := 'diabetes'],
    keys[weight == TRUE][, indication := 'weight_management'],
    keys[!diabetes & !weight][, indication := indication_status]), use.names = TRUE)
}
select_indication_ids <- function(exposures, drug_name, stratum_name) {
  unique(exposures[drug == drug_name & indication == stratum_name, PRIMARYID])
}
run_indication_window_adjusted <- function(baseline_dir = NULL, project_root = 'F:/GLP1-RA/new', output_root = NULL,
                                          prepared_dir = NULL) {
  suppressPackageStartupMessages(library(data.table))
  if (.Platform$OS.type == 'windows') try(suppressWarnings(Sys.setlocale('LC_CTYPE', 'Chinese (Simplified)_China.936')), silent = TRUE)
  old_threads <- getDTthreads(); setDTthreads(4L); on.exit(setDTthreads(old_threads), add = TRUE)
  result_root <- file.path(project_root, '05_\u7ed3\u679c')
  if (is.null(baseline_dir)) {
    candidates <- sort(list.dirs(file.path(result_root, 'new_reanalysis'), recursive = FALSE), decreasing = TRUE)
    complete <- vapply(candidates, function(p) file.exists(file.path(p, 'RUN_STATUS.txt')) &&
      identical(readLines(file.path(p, 'RUN_STATUS.txt'), warn = FALSE)[1], 'COMPLETED'), logical(1))
    if (!any(complete)) stop('No completed baseline')
    baseline_dir <- candidates[complete][1]
  }
  baseline_dir <- normalizePath(baseline_dir, winslash = '/')
  stopifnot(readLines(file.path(baseline_dir, 'RUN_STATUS.txt'), warn = FALSE)[1] == 'COMPLETED')
  sys.source(file.path(baseline_dir, 'code', '01_helpers.R'), envir = environment())
  sys.source(file.path(project_root, '12_firth_reporting.R'), envir = environment())
  config <- readRDS(file.path(baseline_dir, 'analysis_config.rds'))
  if (is.null(output_root)) output_root <- file.path(result_root, 'indication_window_adjusted')
  out <- file.path(output_root, paste0('run_', format(Sys.time(), '%Y%m%d_%H%M%S'), '_', Sys.getpid(), '_firth'))
  for (folder in c('tables', 'audit', 'code', 'models')) dir.create(file.path(out, folder), recursive = TRUE)
  complete <- FALSE; writeLines('RUNNING', file.path(out, 'RUN_STATUS.txt'))
  on.exit(if (!complete) writeLines('FAILED_OR_INTERRUPTED', file.path(out, 'RUN_STATUS.txt')), add = TRUE)
  logfile <- file(file.path(out, 'run_log.txt'), 'wt'); sink(logfile, split = TRUE)
  on.exit({sink(); close(logfile)}, add = TRUE)
  emit <- function(x, name, folder = 'tables') fwrite(x, file.path(out, folder, paste0(name, '.csv')), bom = TRUE, na = '')
  cat('Baseline:', baseline_dir, '\nOutput:', out, '\n'); flush.console()
  plan <- c('Exploratory supplementary analysis; revised after software/numerical validation, not a prospectively registered analysis.',
    'Primary comparisons: semaglutide/tirzepatide vs metformin/SGLT2 in diabetes; orlistat in weight_management. Empagliflozin and other strata are exploratory.',
    'PS exposures and baseline drug mapping; comparator combination products excluded. All 2075 conflicting PRIMARYIDs excluded from all joins.',
    'DRUG and INDI restricted to selected DEMO source quarter; REAC target rederived within that quarter. Missing/unmatched provenance is audited.',
    'INDI linked to the exact PRIMARYID + DRUG_SEQ; dual diabetes/weight records can enter both strata. Known other indications are NOT unknown.',
    'Unknown: no specified indication on ANY eligible sequence for that report/drug; known sequences take precedence. Other-known exposures are audited and excluded from the three requested strata.',
    'Weight dictionary excludes WEIGHT DECREASED, a symptom not unambiguous weight-management intent; dictionary is exported.',
    'For each pair/stratum use inclusive overlap of observed FDA_DATE ranges. This is empirical calendar overlap, not an approved-indication treatment cohort.',
    'Exclude reports with both drugs recorded as PS, irrespective of indication, within each common window. One row per report.',
    'Fixed covariates: age per 10 years, sex, calendar year (linear), reporter HCP vs non-HCP; missing age/sex/reporter explicitly indicated.',
    'Age already harmonized to years upstream. Missing age assigned 50 plus indicator; no automatic variable selection. Remove only exact constant/aliased columns.',
    'Firth logistic for ALL estimable comparisons, no lambda tuning; penalized profile likelihood 95% CI and likelihood-ratio p values.',
    'Complete-case sensitivity repeats crude and adjusted estimates on the SAME complete-case sample. Missing-indicator adjustment does not eliminate missingness bias.',
    'BH correction across all estimable comparisons separately for primary missing-category and complete-case analyses; strata overlap.',
    'Crude ROR uses baseline 0.5 correction for any zero cell; original cell counts preserved. ROR signal: a>=3 and lower95>1.',
    'Zero-comparator-event models are explicitly exploratory and cannot establish comparator safety or equivalence.',
    'Prior ridge outputs were invalidated: outcome-informed lambda revisions are not used for final inference.',
    'Implementation equations: https://search.r-project.org/CRAN/refmans/logistf/html/logistf-package.html')
  writeLines(c(paste('baseline_dir =', baseline_dir), plan), file.path(out, 'ANALYSIS_PLAN.txt'))
  scripts <- c(file.path(project_root, c('10_indication_window_adjusted.R','12_firth_reporting.R','99_indication_window_self_check.R')),
    file.path(baseline_dir, 'code', '01_helpers.R'))
  stopifnot(all(file.copy(scripts, file.path(out, 'code'))))
  emit(data.table(path = scripts, md5 = unname(tools::md5sum(scripts))), 'code_checksums')
  raw_paths <- file.path(config$data_dir, paste0(c('DRUG','INDI','REAC'), '.RData'))
  info <- file.info(raw_paths)
  manifest <- data.table(path = raw_paths, size_bytes = info$size, modified_at = format(info$mtime, tz = 'Asia/Shanghai', usetz = TRUE))
  prior <- fread(file.path(baseline_dir, 'tables/input_manifest.csv'), encoding = 'UTF-8')[match(raw_paths, path)]
  stopifnot(identical(as.numeric(prior$size_bytes), as.numeric(info$size)), identical(prior$modified_at, manifest$modified_at))
  emit(manifest, 'raw_input_manifest')
  if (is.null(prepared_dir)) {
    cat('Reading baseline core:', format(Sys.time()), '\n'); flush.console()
    full_core <- as.data.table(readRDS(file.path(baseline_dir, 'data/analysis_core.rds')))
    conflicts <- fread(file.path(baseline_dir, 'tables/dedup_primaryid_case_conflicts.csv'), colClasses = c(PRIMARYID='character'))
    stopifnot(!anyDuplicated(full_core$PRIMARYID), !anyDuplicated(full_core$CASEID), !anyNA(full_core$TARGET))
    emit(conflicts, 'excluded_conflict_ids')
    population <- data.table(metric = c('baseline_reports','baseline_target_reports','conflicting_ids','conflicting_ids_in_baseline'),
      value = c(nrow(full_core), sum(full_core$TARGET), nrow(conflicts), sum(full_core$PRIMARYID %chin% conflicts$PRIMARYID)))
    core <- full_core[!PRIMARYID %chin% conflicts$PRIMARYID, .(PRIMARYID, CASEID, QUARTER, TARGET, FDA_DATE, AGE_YEARS, SEX_GROUP, REPORTER)]
    rm(full_core); gc(FALSE)
    target <- fread(file.path(baseline_dir, 'data/matched_drug_records.csv'), encoding='UTF-8', colClasses=c(PRIMARYID='character', DRUG_SEQ='character'))
    target <- unique(target[main == TRUE & ROLE_COD == 'PS' & drug %chin% c('semaglutide','tirzepatide'),
      .(PRIMARYID, DRUG_SEQ, DRUGNAME, ROLE_COD, drug)])
    target[, source := 'study_drug']
    hits <- fread(file.path(baseline_dir, 'data/additional_drug_name_matches.csv'), encoding='UTF-8', colClasses=c(PRIMARYID='character',DRUG_SEQ='character'))
    sglt2 <- c('empagliflozin','dapagliflozin','canagliflozin','ertugliflozin','bexagliflozin')
    hits <- hits[ROLE_COD == 'PS' & combination_flag == FALSE & ingredient %chin% c('metformin', sglt2, 'orlistat')]
    single <- hits[ingredient %chin% c('metformin','empagliflozin','orlistat'), .(PRIMARYID,DRUG_SEQ,DRUGNAME,ROLE_COD,drug=ingredient)]
    pooled <- hits[ingredient %chin% sglt2, .(PRIMARYID,DRUG_SEQ,DRUGNAME,ROLE_COD,drug='sglt2_inhibitors')]
    comparators <- unique(rbindlist(list(single,pooled))); comparators[, source := 'active_comparator']
    candidates <- unique(rbindlist(list(target,comparators)))
    candidates[, `:=`(DRUG_SEQ=clean_text(DRUG_SEQ),DRUGNAME=clean_text(DRUGNAME))]
    candidates <- candidates[PRIMARYID %chin% core$PRIMARYID]
    core <- core[PRIMARYID %chin% candidates$PRIMARYID]
    rm(target, hits, single, pooled, comparators); gc(FALSE)
    source_audit <- list()
    raw_selected <- function(name) {
      cat('Reading raw',name,':', format(Sys.time()), '\n'); flush.console()
      raw <- load_one_rdata(file.path(config$data_dir,paste0(name,'.RData')), name)
      if (is.numeric(raw$PRIMARYID)) raw <- raw[PRIMARYID %in% suppressWarnings(as.numeric(core$PRIMARYID))]
      else raw <- raw[PRIMARYID %chin% core$PRIMARYID]
      raw[, PRIMARYID := identifier_text(PRIMARYID)]
      stopifnot('FILE' %in% names(raw))
      raw[, source_quarter := indication_source_quarter(FILE)]
      raw <- merge(raw,core[,.(PRIMARYID,retained_quarter=QUARTER)],by='PRIMARYID')
      source_audit[[name]] <<- raw[, .(table=name, rows=.N, missing_source=sum(is.na(source_quarter)),
        different_quarter=sum(!is.na(source_quarter) & source_quarter != retained_quarter),
        selected_quarter=sum(!is.na(source_quarter) & source_quarter == retained_quarter))]
      excluded <- raw[is.na(source_quarter) | source_quarter != retained_quarter]
      if(nrow(excluded)) emit(excluded,paste0('source_excluded_',name),'audit')
      raw[!is.na(source_quarter) & source_quarter == retained_quarter]
    }
    drug <- raw_selected('DRUG')
    drug[, `:=`(DRUG_SEQ=clean_text(DRUG_SEQ), DRUGNAME=clean_text(DRUGNAME), ROLE_COD=toupper(clean_text(ROLE_COD)))]
    signatures <- unique(drug[ROLE_COD=='PS',.(PRIMARYID,DRUG_SEQ,DRUGNAME,ROLE_COD)])
    candidates_bad <- candidates[!signatures,on=.(PRIMARYID,DRUG_SEQ,DRUGNAME,ROLE_COD)]
    emit(candidates_bad,'baseline_drug_hits_without_selected_source','audit')
    exposures <- merge(candidates,signatures,by=c('PRIMARYID','DRUG_SEQ','DRUGNAME','ROLE_COD'))
    emit(exposures,'verified_drug_exposure_records','audit')
    keys <- unique(exposures[,.(PRIMARYID,DRUG_SEQ,drug,source)])
    emit(keys[is.na(DRUG_SEQ) | !nzchar(DRUG_SEQ)],'missing_drug_sequence','audit')
    keys <- keys[!is.na(DRUG_SEQ) & nzchar(DRUG_SEQ)]
    rm(drug, signatures, candidates, exposures); gc(FALSE)
    indi <- raw_selected('INDI')
    indi[, `:=`(DRUG_SEQ=clean_text(DRUG_SEQ),INDI_PT=toupper(clean_text(INDI_PT)))]
    indi <- indi[keys,on=.(PRIMARYID,DRUG_SEQ),nomatch=0L,allow.cartesian=TRUE]
    indi <- unique(indi[,.(PRIMARYID,DRUG_SEQ,INDI_PT,FILE)])
    emit(indi,'linked_selected_source_indications','audit')
    diabetes <- c('TYPE 2 DIABETES MELLITUS','INSULIN-REQUIRING TYPE 2 DIABETES MELLITUS','DIABETES MELLITUS')
    weight <- c('OBESITY','CENTRAL OBESITY','OVERWEIGHT','WEIGHT CONTROL','WEIGHT MANAGEMENT','WEIGHT REDUCTION','BODY MASS INDEX INCREASED')
    emit(rbindlist(list(data.table(indication='diabetes',term=diabetes),data.table(indication='weight_management',term=weight))), 'indication_dictionary')
    classes <- indication_class_for_key(keys, indi, diabetes, weight)
    # A report/drug is unknown only when all eligible sequences lack a specified indication.
    specified_reports <- unique(classes[indication != 'unknown',.(PRIMARYID,drug)])
    classes <- rbindlist(list(classes[indication != 'unknown'], classes[indication == 'unknown'][!specified_reports,on=.(PRIMARYID,drug)]))
    emit(classes,'drug_specific_indication_exposures','audit')
    emit(indi[,.(records=.N,reports=uniqueN(PRIMARYID)),by=INDI_PT][order(-reports)], 'observed_indication_terms')
    emit(classes[,.(reports=uniqueN(PRIMARYID),drug_sequence_keys=.N),by=.(drug,indication)],'drug_specific_indication_counts')
    emit(keys[,.(reports=uniqueN(PRIMARYID),drug_sequence_keys=.N),by=drug],'drug_exposure_key_counts')
    rm(indi); gc(FALSE)
    reac <- raw_selected('REAC'); reac[, PT := toupper(clean_text(PT))]
    events <- unique(reac[PT==config$target_pt,PRIMARYID])
    core[, baseline_TARGET := TARGET]; core[, TARGET := PRIMARYID %chin% events]
    emit(core[TARGET != baseline_TARGET], 'selected_source_event_changes','audit')
    rm(reac); gc(FALSE)
    core[, `:=`(calendar_year=as.integer(format(FDA_DATE,'%Y')),
      SEX_GROUP=fifelse(SEX_GROUP %chin% c('F','M'),SEX_GROUP,'Unknown'),
      REPORTER_GROUP=fcase(REPORTER %chin% c('MD','PH','HP'),'HCP',REPORTER %chin% c('CN','LW','OT'),'non_HCP',default='Unknown'))]
    core[, complete_covariates := !is.na(AGE_YEARS) & SEX_GROUP!='Unknown' & REPORTER_GROUP!='Unknown']
    emit(rbindlist(source_audit),'source_join_audit')
    emit(population,'analysis_population_summary')
    saveRDS(list(core=core,classes=classes,keys=keys,baseline_dir=baseline_dir,manifest=manifest),file.path(out,'audit','prepared_cohort.rds'),compress=FALSE)
  } else {
    prepared <- readRDS(file.path(prepared_dir,'audit','prepared_cohort.rds'))
    stopifnot(identical(prepared$baseline_dir,baseline_dir), identical(prepared$manifest,manifest))
    core <- prepared$core; classes <- prepared$classes; keys <- prepared$keys
    file.copy(file.path(prepared_dir,'audit','prepared_cohort.rds'), file.path(out,'audit','prepared_cohort.rds'))
    writeLines(prepared_dir,file.path(out,'PREPARED_FROM.txt'))
  }
  emit(core,'analysis_covariate_frame','audit')
  results <- list(); windows <- list(); covariates <- list(); balance <- list()
  for (target_name in c('semaglutide','tirzepatide')) for(comp_name in c('metformin','sglt2_inhibitors','empagliflozin','orlistat')) for(stratum_name in c('diabetes','weight_management','unknown')) {
    label <- paste(target_name,comp_name,stratum_name,sep='__')
    cat('Fitting',label,format(Sys.time()),'\n'); flush.console()
    tpre <- core[PRIMARYID %chin% select_indication_ids(classes,target_name,stratum_name)]
    cpre <- core[PRIMARYID %chin% select_indication_ids(classes,comp_name,stratum_name)]
    tv <- tpre[!is.na(FDA_DATE)]; cv <- cpre[!is.na(FDA_DATE)]
    start <- if(nrow(tv) && nrow(cv)) max(min(tv$FDA_DATE),min(cv$FDA_DATE)) else as.Date(NA)
    end <- if(nrow(tv) && nrow(cv)) min(max(tv$FDA_DATE),max(cv$FDA_DATE)) else as.Date(NA)
    # Exclude any PS coexposure to this pair even if its indication differs.
    both <- intersect(keys[drug==target_name,PRIMARYID],keys[drug==comp_name,PRIMARYID])
    tsel <- tv[!is.na(start) & FDA_DATE>=start & FDA_DATE<=end]
    csel <- cv[!is.na(start) & FDA_DATE>=start & FDA_DATE<=end]
    windows[[label]] <- data.table(comparison=label,target_drug=target_name,comparator=comp_name,indication=stratum_name,
      target_pre=nrow(tpre),comparator_pre=nrow(cpre),target_missing_date=sum(is.na(tpre$FDA_DATE)),comparator_missing_date=sum(is.na(cpre$FDA_DATE)),
      common_start=start,common_end=end,target_in_window=nrow(tsel),comparator_in_window=nrow(csel),
      target_coexposed_excluded=sum(tsel$PRIMARYID %chin% both),comparator_coexposed_excluded=sum(csel$PRIMARYID %chin% both))
    model_data <- rbindlist(list(tsel[!PRIMARYID %chin% both][,EXPOSURE:='target'],csel[!PRIMARYID %chin% both][,EXPOSURE:='comparator']))
    stopifnot(!anyDuplicated(model_data$PRIMARYID))
    covariates[[label]] <- model_data[,.(reports=.N,events=sum(TARGET),age_missing=sum(is.na(AGE_YEARS)),age_missing_pct=100*mean(is.na(AGE_YEARS)),
      sex_missing=sum(SEX_GROUP=='Unknown'),sex_missing_pct=100*mean(SEX_GROUP=='Unknown'),reporter_missing=sum(REPORTER_GROUP=='Unknown'),
      reporter_missing_pct=100*mean(REPORTER_GROUP=='Unknown'),complete_cases=sum(complete_covariates),complete_case_events=sum(TARGET & complete_covariates)),by=EXPOSURE][,comparison:=label]
    if(nrow(model_data)) {
      balance[[label]] <- rbindlist(lapply(c('SEX_GROUP','REPORTER_GROUP','calendar_year'),function(variable) {
        z<-model_data[,.(reports=.N,events=sum(TARGET)),by=c('EXPOSURE',variable)]; setnames(z,variable,'level'); z[,level:=as.character(level)]; z[,variable:=variable]; z
      }))[,comparison:=label]
      age_stats <- model_data[,.(reports=.N,age_known=sum(!is.na(AGE_YEARS)),mean_age=mean(AGE_YEARS,na.rm=TRUE),sd_age=sd(AGE_YEARS,na.rm=TRUE),
        median_age=median(AGE_YEARS,na.rm=TRUE),min_year=min(calendar_year),max_year=max(calendar_year)),by=EXPOSURE]
      emit(age_stats,paste0(label,'__age_summary'),'models')
    }
    for(sample_name in c('missing_category','complete_case')) {
      dd <- if(sample_name=='complete_case') model_data[complete_covariates==TRUE] else model_data
      a<-sum(dd$TARGET & dd$EXPOSURE=='target'); b<-sum(!dd$TARGET & dd$EXPOSURE=='target')
      cc<-sum(dd$TARGET & dd$EXPOSURE=='comparator'); d<-sum(!dd$TARGET & dd$EXPOSURE=='comparator')
      row <- signal_metrics(a,b,cc,d)
      row[,`:=`(comparison=label,target_drug=target_name,comparator=comp_name,indication=stratum_name,sample=sample_name,common_start=start,common_end=end,
        primary_comparison=(stratum_name=='diabetes' && comp_name %chin% c('metformin','sglt2_inhibitors')) || (stratum_name=='weight_management' && comp_name=='orlistat'))]
      if(a+b>0 && cc+d>0 && a+cc>0 && b+d>0) {
        fit <- fit_adjusted_ror(dd)
        emit(fit$grouped,paste0(label,'__',sample_name,'__grouped'),'models')
        emit(fit$coefficients,paste0(label,'__',sample_name,'__coefficients'),'models')
        row <- cbind(row,fit$summary)
      } else row[,`:=`(adjusted_model='not_fit',adjusted_model_status='no_two_groups_or_no_outcome_variation',adjusted_n=nrow(dd),adjusted_events=a+cc,
        adjusted_ROR=NA_real_,adjusted_ROR_low=NA_real_,adjusted_ROR_high=NA_real_,adjusted_p=NA_real_)]
      row[,`:=`(zero_comparator_events=cc==0, sparse_either_group=min(a,cc)<10,
        adjusted_signal=ifelse(is.na(adjusted_ROR_low),NA,a>=3 & adjusted_ROR_low>1))]
      results[[paste(label,sample_name)]] <- row
    }
  }
  results <- rbindlist(results,fill=TRUE)
  results[, `:=`(crude_ROR=ROR,crude_ROR_low=ROR_low,crude_ROR_high=ROR_high,crude_ROR_signal=ROR_signal)]
  results[,`:=`(adjusted_q_BH=p.adjust(adjusted_p,'BH'),fisher_q_BH=p.adjust(fisher_p,'BH')),by=sample]
  stopifnot(all(results$total_reports==results$a+results$b+results$c+results$d),
    all(results$total_reports==results$adjusted_n),all(results$adjusted_events==results$a+results$c),
    all(results$profile_endpoint_error[!is.na(results$profile_endpoint_error)]<1e-3))
  emit(results,'indication_common_window_adjusted_results')
  emit(rbindlist(windows),'common_calendar_window_audit')
  emit(rbindlist(covariates,fill=TRUE),'covariate_missingness_by_group')
  emit(rbindlist(balance,fill=TRUE),'covariate_distributions_by_group')
  emit(results[,.(comparison,sample,adjusted_n,adjusted_events,parameters,events_per_parameter,adjusted_model,adjusted_model_status,
    max_score,iterations,profile_endpoint_error,removed_constant_or_alias_columns)],'adjusted_model_audit')
  capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
  writeLines('COMPLETED_AWAITING_INDEPENDENT_VALIDATION',file.path(out,'RUN_STATUS.txt')); complete<-TRUE
  cat('MODELS COMPLETED:',out,'\n'); invisible(out)
}
if(sys.nframe()==0L) indication_window_output <- run_indication_window_adjusted()