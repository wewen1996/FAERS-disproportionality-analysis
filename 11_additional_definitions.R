additional_definitions <- function() {
  make_aliases <- function(ingredient, category, subclass, aliases) {
    data.table::data.table(ingredient = ingredient, category = category,
      subclass = subclass, alias = unique(toupper(aliases)))
  }
  aliases <- list(
    make_aliases("metformin", "glucose_lowering", "biguanide", c("METFORMIN", "GLUCOPHAGE", "GLUMETZA", "FORTAMET", "RIOMET")),
    make_aliases("empagliflozin", "glucose_lowering", "SGLT2", c("EMPAGLIFLOZIN", "JARDIANCE")),
    make_aliases("dapagliflozin", "glucose_lowering", "SGLT2", c("DAPAGLIFLOZIN", "FARXIGA", "FORXIGA")),
    make_aliases("canagliflozin", "glucose_lowering", "SGLT2", c("CANAGLIFLOZIN", "INVOKANA")),
    make_aliases("ertugliflozin", "glucose_lowering", "SGLT2", c("ERTUGLIFLOZIN", "STEGLATRO")),
    make_aliases("bexagliflozin", "glucose_lowering", "SGLT2", c("BEXAGLIFLOZIN", "BRENZAVVY")),
    make_aliases("orlistat", "weight_management", "lipase_inhibitor", c("ORLISTAT", "XENICAL", "ALLI")),
    make_aliases("furosemide", "diuretic", "loop", c("FUROSEMIDE", "FRUSEMIDE", "LASIX")),
    make_aliases("bumetanide", "diuretic", "loop", c("BUMETANIDE", "BUMEX")),
    make_aliases("torsemide", "diuretic", "loop", c("TORSEMIDE", "TORASEMIDE", "DEMADEX")),
    make_aliases("hydrochlorothiazide", "diuretic", "thiazide", c("HYDROCHLOROTHIAZIDE", "HCTZ", "MICROZIDE")),
    make_aliases("spironolactone", "diuretic", "potassium_sparing", c("SPIRONOLACTONE", "ALDACTONE")),
    make_aliases("eplerenone", "diuretic", "potassium_sparing", c("EPLERENONE", "INSPRA")),
    make_aliases("warfarin", "anticoagulant", "VKA", c("WARFARIN", "COUMADIN", "JANTOVEN")),
    make_aliases("apixaban", "anticoagulant", "DOAC", c("APIXABAN", "ELIQUIS")),
    make_aliases("rivaroxaban", "anticoagulant", "DOAC", c("RIVAROXABAN", "XARELTO")),
    make_aliases("dabigatran", "anticoagulant", "DOAC", c("DABIGATRAN", "PRADAXA")),
    make_aliases("edoxaban", "anticoagulant", "DOAC", c("EDOXABAN", "SAVAYSA", "LIXIANA")),
    make_aliases("enoxaparin", "anticoagulant", "heparin_related", c("ENOXAPARIN", "LOVENOX", "CLEXANE")),
    make_aliases("dalteparin", "anticoagulant", "heparin_related", c("DALTEPARIN", "FRAGMIN")),
    make_aliases("insulin", "glucose_lowering", "insulin", c("INSULIN", "HUMALOG", "NOVOLOG", "NOVORAPID", "LANTUS", "LEVEMIR", "TRESIBA", "HUMULIN", "NOVOLIN", "TOUJEO", "BASAGLAR", "FIASP", "APIDRA"))
  )
  generic_sets <- list(
    diuretic = c("CHLOROTHIAZIDE", "CHLORTHALIDONE", "INDAPAMIDE", "METOLAZONE", "BENDROFLUMETHIAZIDE",
      "AMILORIDE", "TRIAMTERENE", "ETHACRYNIC ACID", "AZOSEMIDE", "ACETAZOLAMIDE", "MANNITOL"),
    anticoagulant = c("HEPARIN", "TINZAPARIN", "NADROPARIN", "PARNAPARIN", "BEMIPARIN", "FONDAPARINUX",
      "ACENOCOUMAROL", "PHENPROCOUMON", "FLUINDIONE", "ARGATROBAN", "BIVALIRUDIN", "DANAPAROID"),
    glucose_lowering = c("GLIPIZIDE", "GLIMEPIRIDE", "GLICLAZIDE", "GLIBENCLAMIDE", "GLYBURIDE",
      "TOLBUTAMIDE", "CHLORPROPAMIDE", "SITAGLIPTIN", "LINAGLIPTIN", "SAXAGLIPTIN", "ALOGLIPTIN",
      "VILDAGLIPTIN", "TENELIGLIPTIN", "PIOGLITAZONE", "ROSIGLITAZONE", "ACARBOSE", "MIGLITOL",
      "VOGLIBOSE", "REPAGLINIDE", "NATEGLINIDE"))
  for (category in names(generic_sets)) {
    for (ingredient in generic_sets[[category]]) {
      aliases[[length(aliases) + 1L]] <- make_aliases(tolower(gsub(" ", "_", ingredient)), category, "other", ingredient)
    }
  }
  aliases <- data.table::rbindlist(aliases)
  sglt2 <- c("empagliflozin", "dapagliflozin", "canagliflozin", "ertugliflozin", "bexagliflozin")
  comparators <- list(metformin = "metformin", empagliflozin = "empagliflozin",
    sglt2_inhibitors = sglt2, orlistat = "orlistat")
  indications <- list(
    type2_diabetes = c("TYPE 2 DIABETES MELLITUS", "INSULIN-REQUIRING TYPE 2 DIABETES MELLITUS"),
    diabetes_broad = c("TYPE 2 DIABETES MELLITUS", "INSULIN-REQUIRING TYPE 2 DIABETES MELLITUS", "DIABETES MELLITUS"),
    weight_management = c("OBESITY", "CENTRAL OBESITY", "OVERWEIGHT", "WEIGHT CONTROL", "WEIGHT MANAGEMENT",
      "WEIGHT REDUCTION", "BODY MASS INDEX INCREASED", "WEIGHT DECREASED"))
  conditions <- list(
    vascular_disease = c("CORONARY ARTERY DISEASE", "CORONARY ARTERIOSCLEROSIS", "MYOCARDIAL ISCHAEMIA",
      "ISCHAEMIC HEART DISEASE", "ANGINA PECTORIS", "MYOCARDIAL INFARCTION", "ATHEROSCLEROSIS",
      "ARTERIOSCLEROSIS", "PERIPHERAL ARTERIAL OCCLUSIVE DISEASE", "PERIPHERAL ARTERY DISEASE",
      "PERIPHERAL VASCULAR DISORDER", "CEREBROVASCULAR ACCIDENT", "ISCHAEMIC STROKE", "TRANSIENT ISCHAEMIC ATTACK"),
    heart_failure = c("CARDIAC FAILURE", "CARDIAC FAILURE CONGESTIVE", "CARDIAC FAILURE CHRONIC",
      "CARDIAC FAILURE ACUTE", "LEFT VENTRICULAR FAILURE", "RIGHT VENTRICULAR FAILURE", "CONGESTIVE HEART FAILURE",
      "HEART FAILURE", "HEART FAILURE WITH REDUCED EJECTION FRACTION", "HEART FAILURE WITH PRESERVED EJECTION FRACTION",
      "HEART FAILURE WITH MIDRANGE EJECTION FRACTION", "HEART FAILURE WITH IMPROVED EJECTION FRACTION"),
    constipation = c("CONSTIPATION", "CHRONIC CONSTIPATION", "FAECAL IMPACTION", "FAECALOMA"),
    dehydration_hypovolaemia = c("DEHYDRATION", "HYPOVOLAEMIA"),
    hypotension = c("HYPOTENSION", "ORTHOSTATIC HYPOTENSION"),
    diabetes_recorded = indications$diabetes_broad,
    weight_management_recorded = indications$weight_management,
    metabolic_or_weight_recorded = unique(c(indications$diabetes_broad, indications$weight_management)))
  dates <- data.table::data.table(
    node_id = c("reference_temporal_proxy", "sema_weight_approval", "tirz_weight_approval", "sema_colitis_case", "tirz_colitis_case"),
    drug = c(NA_character_, "semaglutide", "tirzepatide", "semaglutide", "tirzepatide"),
    date = as.Date(c("2023-07-01", "2021-06-04", "2023-11-08", "2025-02-10", "2024-11-06")),
    node_type = c("reference_paper_temporal_proxy", "primary_regulatory_attention_proxy", "primary_regulatory_attention_proxy",
      "secondary_colitis_publication_proxy", "secondary_colitis_publication_proxy"),
    source_id = c("REFERENCE_PAPER_TEMPORAL", "FDA_WEGOVY_2021", "FDA_ZEPBOUND_2023", "SEMA_CASE_2025", "TIRZ_CASE_2024"))
  # Keep each citation in one record; unequal column vectors can recycle and
  # silently associate a source_id with the next source's identifier or URL.
  source_record <- function(source_id, title, identifier, url) {
    fields <- list(source_id = source_id, title = title, identifier = identifier, url = url)
    stopifnot(all(lengths(fields) == 1L), !anyNA(unlist(fields)), all(nzchar(unlist(fields))))
    data.table::as.data.table(c(fields, list(checked_on = "2026-10-13")))
  }
  sources <- data.table::rbindlist(list(
    source_record("REFERENCE_PAPER",
      "Exploring potential associations between GLP-1RAs and depressive disorders: Table 5 and temporal methods",
      "10.1016/j.eclinm.2025.103385", "https://doi.org/10.1016/j.eclinm.2025.103385"),
    source_record("REFERENCE_PAPER_TEMPORAL",
      "Reference paper temporal stratification: Q2 2023 or earlier versus Q3 2023 onwards; exploratory proxy, not an ischemic-colitis warning date",
      "10.1016/j.eclinm.2025.103385", "https://doi.org/10.1016/j.eclinm.2025.103385"),
    source_record("FDA_WEGOVY_2021", "Wegovy original FDA approval letter", "NDA215256",
      "https://www.accessdata.fda.gov/drugsatfda_docs/appletter/2021/215256Orig1s000ltr.pdf"),
    source_record("FDA_ZEPBOUND_2023", "FDA approves Zepbound for chronic weight management", "NDA217806",
      "https://www.fda.gov/news-events/press-announcements/fda-approves-new-medication-chronic-weight-management"),
    source_record("FDA_OZEMPIC_2017", "Drugs@FDA NDA209637 original submission", "NDA209637",
      "https://api.fda.gov/drug/drugsfda.json?search=application_number:NDA209637&limit=1"),
    source_record("FDA_MOUNJARO_2022", "Drugs@FDA NDA215866 original submission", "NDA215866",
      "https://api.fda.gov/drug/drugsfda.json?search=application_number:NDA215866&limit=1"),
    source_record("SEMA_CASE_2025",
      "Semaglutide-Associated Ischemic Colitis in a Patient Without Traditional Risk Factors: A Case Report",
      "10.7759/cureus.78832; PMID 40084311; PMC11904261", "https://pmc.ncbi.nlm.nih.gov/articles/PMC11904261/"),
    source_record("TIRZ_CASE_2024", "Tirzepatide-Associated Colonic Ischemia",
      "10.14309/crj.0000000000001551; PMID 39507503; PMC11540430", "https://pmc.ncbi.nlm.nih.gov/articles/PMC11540430/")
  ))
  stopifnot(!anyDuplicated(sources$source_id), all(dates$source_id %in% sources$source_id))
  list(aliases = aliases, comparators = comparators, indications = indications, conditions = conditions,
    nodes = dates, sources = sources,
    first_us_approval = as.Date(c(semaglutide = "2017-12-05", tirzepatide = "2022-05-13")),
    combination_brand_pattern = "(^| )(SYNJARDY|XIGDUO|INVOKAMET|SEGLUROMET|TRIJARDY|GLYXAMBI|QTERN|STEGLUJAN|JANUMET|JENTADUETO|KOMBIGLYZE|KAZANO|EUCREAS|GALVUS MET|ACTOPLUS MET|AVANDAMET|GLUCOVANCE|METAGLIP)( |$)",
    hill_criteria = c("Strength of the association", "Consistency", "Analogy", "Biological plausibility/empirical evidence",
      "Exclusion of biases/confounders", "Specificity", "Temporal relationship", "Reversibility", "Coherence"))
}

normalize_drug_words <- function(value) {
  trimws(gsub(" +", " ", gsub("[^A-Z0-9]+", " ", toupper(clean_text(value)))))
}

match_additional_drug_names <- function(catalog, definitions) {
  catalog <- data.table::copy(catalog)
  catalog[, `:=`(NAME_WORDS = normalize_drug_words(DRUGNAME), AI_WORDS = normalize_drug_words(PROD_AI))]
  forbidden <- grepl("(^| )(PLACEBO|BLINDED|DUMMY)( |$)", catalog$NAME_WORDS)
  matched <- list()
  for (ingredient_name in unique(definitions$aliases$ingredient)) {
    aliases <- definitions$aliases[ingredient == ingredient_name]
    pattern <- paste0("(^| )(", paste(normalize_drug_words(aliases$alias), collapse = "|"), ")( |$)")
    match_name <- grepl(pattern, catalog$NAME_WORDS)
    match_ai <- grepl(pattern, catalog$AI_WORDS)
    selected <- which((match_name | match_ai) & !forbidden)
    if (length(selected)) matched[[ingredient_name]] <- data.table::data.table(
      NAME_ID = catalog$NAME_ID[selected], ingredient = ingredient_name,
      category = aliases$category[1L], subclass = aliases$subclass[1L],
      match_source = ifelse(match_name[selected] & match_ai[selected], "both", ifelse(match_ai[selected], "PROD_AI", "DRUGNAME")))
  }
  hits <- if (length(matched)) data.table::rbindlist(matched) else data.table::data.table(
    NAME_ID = integer(), ingredient = character(), category = character(), subclass = character(), match_source = character())
  catalog[, combination_flag := grepl(";|/|\\+| AND ", toupper(clean_text(PROD_AI))) |
    grepl(definitions$combination_brand_pattern, NAME_WORDS)]
  multiple <- hits[, .(ingredient_n = data.table::uniqueN(ingredient)), by = NAME_ID][ingredient_n > 1L, NAME_ID]
  catalog[NAME_ID %in% multiple, combination_flag := TRUE]
  hits <- merge(hits, catalog[, .(NAME_ID, combination_flag)], by = "NAME_ID", all.x = TRUE)
  list(hits = hits, catalog = catalog)
}

proportion_recorded <- function(n, denominator) {
  stopifnot(n >= 0, denominator >= n)
  if (!denominator) return(c(percent = NA_real_, low = NA_real_, high = NA_real_))
  probability <- n / denominator
  z_value <- stats::qnorm(0.975)
  correction <- 1 + z_value^2 / denominator
  midpoint <- (probability + z_value^2 / (2 * denominator)) / correction
  half_width <- z_value * sqrt(probability * (1 - probability) / denominator + z_value^2 / (4 * denominator^2)) / correction
  c(percent = 100 * probability, low = 100 * max(0, midpoint - half_width), high = 100 * min(1, midpoint + half_width))
}
