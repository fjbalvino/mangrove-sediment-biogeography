# Shared, frozen broad-category keyword rules.
# Used unchanged for primary dictionary-only annotation and partial-MODULE sensitivity.
FUNCTIONAL_CATEGORY_ORDER <- c(
  "Anaerobic redox / energy",
  "Transport / membrane systems",
  "Information processing",
  "Amino acid / cofactor biosynthesis",
  "Carbon / carbohydrate metabolism",
  "Signal / stress / metals",
  "Uncharacterized / other"
)

FUNCTIONAL_PATTERNS <- list(
  "Anaerobic redox / energy" = c(
    "ferredoxin", "oxidoreductase", "reductase", "dehydrogenase",
    "hydrogenase", "heterodisulfide", "formate", "carbon-monoxide",
    "carbon monoxide", "co dehydrogenase", "codh", "acetyl-coa synthase",
    "qmo", "hdr", "dms", "dimethyl sulfoxide", "nitrate", "nitrite",
    "sulfate", "sulfite", "sulfur", "methan", "nadh-quinone",
    "succinate dehydrogenase", "fumarate reductase", "electron transfer"
  ),
  "Transport / membrane systems" = c(
    "transport", "abc", "permease", "substrate-binding", "binding protein",
    "efflux", "membrane fusion", "transporter", "translocation",
    "nickel transport", "peptide transport", "sugar transport"
  ),
  "Information processing" = c(
    "trna", "synthetase", "ribosomal", "ribosome", "translation",
    "transcription", "polymerase", "elongation factor", "nusa",
    "uvr", "mutl", "muts", "repair", "replication", "recombin",
    "preprotein translocase", "seca"
  ),
  "Amino acid / cofactor biosynthesis" = c(
    "leucine", "isoleucine", "valine", "aspartate", "homoserine",
    "threonine", "lysine", "dap", "amino acid", "aminotransferase",
    "molybdopterin", "nad", "cofactor", "pantothenate", "biotin",
    "folate", "methylenetetrahydrofolate", "carbamoyl"
  ),
  "Carbon / carbohydrate metabolism" = c(
    "glycolysis", "gluconeogenesis", "sugar", "ribose", "galactose",
    "carbohydrate", "acetate", "pyruvate", "lactate", "acetoin",
    "coa ligase", "phosphofructokinase", "glucan", "starch",
    "cellulose", "glycosidase"
  ),
  "Signal / stress / metals" = c(
    "two-component", "response regulator", "sensor kinase", "copper",
    "silver", "metal", "selenocysteine", "stress", "chaperone",
    "heat shock", "cold shock", "detox", "resistance"
  )
)

score_functional_category <- function(function_text, module_text) {
  text <- tolower(paste(
    ifelse(is.na(function_text), "", function_text),
    ifelse(is.na(module_text), "", module_text)
  ))
  scores <- vapply(FUNCTIONAL_PATTERNS, function(patterns) {
    sum(vapply(patterns, function(pattern) {
      grepl(pattern, text, fixed = TRUE)
    }, logical(1)))
  }, numeric(1))
  max_score <- max(scores)
  if (max_score == 0L) {
    selected <- "Uncharacterized / other"
    tied <- FALSE
    candidates <- NA_character_
  } else {
    top <- names(scores)[scores == max_score]
    # FUNCTIONAL_PATTERNS order is the documented deterministic tie breaker.
    selected <- top[1L]
    tied <- length(top) > 1L
    candidates <- paste(top, collapse = "; ")
  }
  list(selected = selected, max_score = max_score, tied = tied,
       candidates = candidates, scores = scores)
}

