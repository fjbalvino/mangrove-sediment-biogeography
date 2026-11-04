"""Use the shared hierarchical summaries and descriptive models."""
from pathlib import Path
import importlib.util
p=Path(__file__).resolve().parents[2]/'analysis_nitrogen/scripts/12_gene_environment/12B_BUILD_NITROGEN_PILOT.py'
s=importlib.util.spec_from_file_location('shared_gene_models',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
load_inputs=m.load_inputs
hierarchical_summary=m.hierarchical_summary
association_diagnostics=m.association_diagnostics
