# Publishing the interactive atlas

## Recommended architecture

Keep the analysis and viewer generators in this manuscript repository. Generate and inspect the complete website locally, then deploy a reviewed website artifact through **GitHub Pages + GitHub Actions**. Do not commit generated viewers, matrices, images or simulation outputs to a branch, including a `gh-pages` branch.

The website artifact can be stored as a GitHub Release attachment uploaded from the local computer; it does not enter Git history. This still uploads content to GitHub: the HTML embeds selected derived gene values, coordinates, environmental estimates and rendered images. Publishing a functioning atlas necessarily exposes these display values. Encoding them in HTML is not a way to call them “code only”. Under the current instruction, this update publishes only the source and this procedure, **not the generated website or Release asset**.

GitHub supports ordinary images; keeping them out of this repository is a project organization/reproducibility choice, not a blanket GitHub prohibition. Large file/tool transport failures should not be solved by putting encoded matrices or images into source files.

| Route | Suitability |
| --- | --- |
| Local atlas, linked from repository documentation | Current choice until publication is authorized; complete interactivity available locally |
| Pages deployed from a locally reviewed Release artifact | Recommended public route: website remains separate from source history; repeatable version/hash verification |
| Commit the complete archive or generated HTML to Git | Avoid: mixes data/images with code and duplicates large embedded payloads across history |
| External institutional web hosting plus a repository link | Appropriate if the website should not be uploaded to GitHub at all; the same static export can be used |

## Export locally

After generating the run's atlas:

```bash
python workflow/09_EXPORT_ATLAS_SITE.py \
  --atlas /path/to/runs/functional-complete/atlas \
  --out /path/to/publication/atlas-site
```

The exporter creates `atlas-site/` and `atlas-site.zip`, checks internal file links and rejects symlinks, unexpected file types and a site reaching 1 GB. It copies the index, eight self-contained viewers and licensing/manifest text. It excludes complete metadata, raw KO matrices, numerical result tables, checkpoints and manuscript PDF/Word files. Embedded display data and downloadable heatmap images/sample-order information remain part of the viewers.

Review the actual website before any upload:

```bash
python -m http.server 8000 --directory /path/to/publication/atlas-site
```

Open `http://localhost:8000`, test all eight panels, zoom to overlapping samples, change genes/thresholds, inspect heatmaps and environmental sensitivity, and disable terrain to inspect the offline fallback. External map tiles require an internet connection. The local HTTP server is only a preview, not public hosting.

## Publish later from the local computer

The following steps deliberately upload the reviewed site and should be used only when its publication is authorized. They have **not** been run by this update.

1. Copy `workflow/atlas_pages.yml.example` to `.github/workflows/atlas-pages.yml`, commit the workflow text, and in repository **Settings → Pages → Source** choose **GitHub Actions**. The example is inactive in its current location.
2. Authenticate locally with GitHub CLI (`gh auth login`). Create a dedicated, descriptive release tag and upload only the site ZIP as a Release asset. For example, replace `atlas-v1` if that tag is already used:

```bash
gh release create atlas-v1 /path/to/publication/atlas-site.zip \
  --repo fjbalvino/mangrove-sediment-biogeography --target main \
  --title "Biogeochemical gene atlas" \
  --notes "Reviewed static atlas. See the included site manifest and repository reproduction guide."
sha256sum /path/to/publication/atlas-site.zip
```

3. Trigger the workflow from the local terminal, substituting the exact 64-character digest printed above:

```bash
gh workflow run atlas-pages.yml \
  --repo fjbalvino/mangrove-sediment-biogeography --ref main \
  -f release_tag=atlas-v1 -f sha256=REPLACE_WITH_EXACT_SHA256
```

4. Review the Actions run and the deployment URL it reports. Only after confirming it works, add that URL to the repository README and manuscript atlas reference. A likely project URL is not evidence that a site is deployed.

The workflow downloads only the named Release asset, checks its SHA-256 and internal paths, validates the site, and uses the official Pages upload/deploy actions. It does not execute analysis code from the ZIP or import the full evidence archive. Generated website content stays outside the repository's Git history. Static content is public after release/deployment.

## Current platform documentation

Checked 3 November 2026:

- [GitHub Pages custom workflows](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages): artifact-based deployment and required permissions/environment. The example follows `configure-pages@v5`, `upload-pages-artifact@v4`, `deploy-pages@v4` and `checkout@v6` shown in the documentation.
- [GitHub Pages limits](https://docs.github.com/en/pages/getting-started-with-github-pages/github-pages-limits): published site at most 1 GB; soft monthly bandwidth limit 100 GB. The published-site limit applies even though the artifact transport supports larger archives.
- [GitHub CLI release create](https://cli.github.com/manual/gh_release_create) and [workflow run](https://cli.github.com/manual/gh_workflow_run): local upload and manual deployment commands.

For a substantially larger future dataset, reassess payload size and loading time before publishing. Separate shared visual assets and per-cycle derived data if needed; maintain the same scientific quantities and missing-data semantics. Do not change the analysis merely to fit a hosting limit.
