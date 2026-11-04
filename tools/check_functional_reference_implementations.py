#!/usr/bin/env python3
"""Execute small probes of pinned upstream functions, never study data.

Supply local checkouts at SOURCES/owner/repository. No automatic download and no
third-party source redistribution. This is NOT an end-to-end package benchmark.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    if not __debug__:
        raise RuntimeError('Do not run with Python -O; probe assertions must remain active.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    sources = args.sources.resolve()
    manifest = json.loads((ROOT/'config/gene_biogeography/upstream_sources.json').read_text())
    lookup = {(r['repository'], f['path']): f for r in manifest['repositories'] for f in r['files']}
    used = []

    def source(repo, path):
        p = sources/repo/path
        data = p.read_bytes()
        sha = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        if sha != lookup[repo,path]['git_blob_sha']:
            raise ValueError('Upstream source mismatch: '+str(p))
        used.append(dict(repository=repo,path=path,git_blob_sha=sha))
        return p

    results = []
    with tempfile.TemporaryDirectory() as temp:
        temp = Path(temp)
        # The upstream Pathway class has no package/runtime dependencies.
        path = source('SilentGene/DiTing','diting/workflow/scripts/pathway.py')
        diting = load_module('diting_probe',path)
        obj = diting.Pathway({'S_K02588':9.0}, 'S')
        obj.Nitrogen_fixation()
        val = next(iter(obj.out_data.values()))
        assert val == 3 and obj.abundance['S_K02586'] == 0
        results.append(dict(probe='DiTing_partial_nitrogenase',value=val,
            finding='One subunit produces a positive arithmetic score; missing dictionary entries become zero.'))

        # Run the complete upstream filter CLI, with a below-model-threshold hit.
        filt = source('SilentGene/DiTing','diting/workflow/scripts/kofamscan_filter.py')
        input_path = temp/'hits.tsv'; output_path = temp/'hits_filtered.tsv'
        input_path.write_text('gene1\tK02588\t100\t50\t1e-20\tnitrogenase\n')
        subprocess.run([sys.executable,str(filt),'-i',str(input_path),'-o',str(output_path),'-s','S'],check=True,capture_output=True,text=True)
        result = output_path.read_text().strip().split('\t')
        assert result == ['S','gene1','K02588']
        results.append(dict(probe='DiTing_subthreshold_filter',output_columns=len(result),
            finding='Below-threshold hit retained; exported row omits score and threshold flag.'))

        # Run the complete PCycDB coverage CLI on deliberately overlapping names.
        cov_script = source('ZengJiaxiong/Phosphorus-cycling-database','Coverage_get.py')
        cov = temp/'coverage.tsv'; mapping = temp/'mapping.tsv'; out = temp/'coverage_out.tsv'
        cov.write_text('orf1\tx\tx\t2\norf2\tx\tx\t3\n')
        mapping.write_text('orf1\tfoo\norf2\tfooB\n')
        subprocess.run([sys.executable,str(cov_script),'-i',str(mapping),'-t',str(cov),'-o',str(out)],check=True,capture_output=True,text=True)
        observed = dict(line.split('\t') for line in out.read_text().splitlines())
        assert float(observed['foo']) == 5 and float(observed['fooB']) == 3
        results.append(dict(probe='PCycDB_prefix_collision',observed=observed,
            finding='startswith aggregates fooB into foo; exact identifiers are required in any adapter.'))

        # Unmodified DRAM score helper; this does not execute annotation searches.
        dram = load_module('dram_hmm_probe',source('BortonWrightonLabs/DRAM','bin/hmm_parser.py'))
        import pandas as pd
        frame = pd.DataFrame({'score_type':['full','domain'],'full_bitScore':[9.,9.],
                              'domain_bitScore':[100.,12.],'threshold':[10.,10.]})
        assert dram.sig_scores_row_by_row(frame,'kofam').index.tolist() == [1]
        results.append(dict(probe='DRAM_score_type',finding='Full/domain score choice and threshold rejection verified.'))

        # Extract unmodified R function expressions, avoiding package startup.
        cnps = source('yuezhengfu/CNPS.cycle','R/Ncycleabundance.R')
        maaslin = source('biobakery/maaslin3','R/utility_scripts.R')
        ggkegg = source('noriakis/ggkegg','R/module_functions.R')
        rscript = temp/'probe.R'
        rscript.write_text(r'''
args <- commandArgs(TRUE)
select_function <- function(file, name) {
  expressions <- parse(file)
  for (e in expressions)
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
        identical(as.character(e[[2]]), name)) {
      eval(e, envir=.GlobalEnv)
      return(invisible(NULL))
    }
  stop(paste("Function not found:",name))
}
select_function(args[1], "Ncyc.abundance")
input <- matrix(9, nrow=1, ncol=2, dimnames=list("K02588",c("S1","S2")))
err <- tryCatch({Ncyc.abundance(input);FALSE},error=function(e) TRUE)
stopifnot(err)
ko <- input
ignored <- input * 10
v <- Ncyc.abundance(ignored)
stopifnot(as.numeric(v["NF.N2.NH4","S1"]) == 3)
cat("CNPS_global_argument_and_partial_score PASS\n")
select_function(args[2], "TSSnorm")
x <- data.frame(A=c(2,0),B=c(2,4),C=c(6,0))
norm <- TSSnorm(x,0)
stopifnot(norm$A[1]==0.2, norm$B[1]==0.2, norm$C[1]==0.6,
          is.na(norm$A[2]), norm$B[2]==1, is.na(norm$C[2]))
cat("MaAsLin_TSS_positive_values_and_mask PASS\n")
select_function(args[3], "module_completeness")
# Only the output constructor is adapted to base data.frame.
# The upstream evaluation function body remains unchanged.
tibble <- function(...) data.frame(...,check.names=FALSE)
setClass("AuditModule",slots=c(definitions="list"))
exprs <- c("K00001+K00002","K00001,K00002",
           "(K00001+K00002),K00003","K00001-K00002")
results <- list()
for (expr in exprs) for (bits in 0:7) {
  allko <- c("K00001","K00002","K00003")
  mask <- as.logical(intToBits(bits)[1:3])
  in_expr <- regmatches(expr,gregexpr("K[0-9]{5}",expr))[[1]]
  block <- list(definition_ko_in_block=list(in_expr),
                definition_block=expr,definition_num_in_block=length(in_expr))
  kmo <- new("AuditModule",definitions=list("1"=block))
  value <- module_completeness(kmo,allko[mask])$complete
  results[[length(results)+1]] <- data.frame(expression=expr,bits=bits,value=as.integer(value))
}
write.table(do.call(rbind,results),args[4],sep="\t",row.names=FALSE,quote=FALSE)
cat("ggkegg_binary_logic_32_cases PASS\n")
''')
        truth = temp/'truth.tsv'
        proc = subprocess.run(['Rscript','--vanilla',str(rscript),str(cnps),str(maaslin),str(ggkegg),str(truth)],
                              check=True,capture_output=True,text=True)
        import numpy as np
        own = load_module('own_rule_probe',ROOT/'scripts/12_gene_environment/template/evaluate.py')
        for row in pd.read_csv(truth,sep='\t').itertuples():
            values = {f'K0000{i+1}': np.array([(row.bits >> i)&1]) for i in range(3)}
            assert own.evaluate(own.parse(row.expression),values,1)[0] == row.value
        results.extend([
            dict(probe='CNPS_global_argument',finding='Without global ko the supplied argument fails; with global ko the supplied argument is ignored; nifH alone scores 3.'),
            dict(probe='MaAsLin_TSS_helper',finding='Matches full-matrix TSS for positives; zeros masked to NA in upstream abundance helper. Not a full MaAsLin fit.'),
            dict(probe='ggkegg_binary_rule_reference',cases=32,finding='All known-row binary cases agree; base data.frame used only for returned table. Unknown-row behavior intentionally differs.')])

        # Evaluate the original Escher constants/helpers in a Node VM.
        escher = source('opencobra/escher','src/dataStyles.js')
        lines = escher.read_text().splitlines()
        snippet = '\n'.join(lines[4:26]+lines[349:450])
        js = temp/'probe.js'
        js.write_text(snippet+'\n'+r'''
const v = {A:[9],B:[null]};
const result = {
  mean:evaluate_gene_reaction_rule("A and B",v,"mean")[0],
  min:evaluate_gene_reaction_rule("A and B",v,"min")[0],
  allMissing:evaluate_gene_reaction_rule("A and B",{A:[null],B:[null]},"min")[0],
  alternative:evaluate_gene_reaction_rule("A or B",{A:[9],B:[3]},"min")[0]
};
if (result.mean !== 4.5 || result.min !== 0 || result.allMissing !== null || result.alternative !== 12) throw Error("Unexpected reference behavior");
process.stdout.write(JSON.stringify(result));
''')
        proc = subprocess.run(['node',str(js)],check=True,capture_output=True,text=True)
        results.append(dict(probe='Escher_GPR_display',observed=json.loads(proc.stdout),
            finding='Partly missing AND input becomes zero before min/mean; all-missing stays null. Display score is not biochemical flux.'))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps({'status':'PASS','scope':'Pinned helper/CLI probes on synthetic inputs only',
                                      'probes':results,'verified_source_files':used},indent=2)+'\n')
    print('REFERENCE_IMPLEMENTATION_PROBES_PASS:',len(results))


if __name__ == '__main__':
    main()
