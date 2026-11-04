#!/usr/bin/env python3
"""Export only the static atlas website for local review and later deployment.

HTML embeds selected derived gene values, coordinates and figure images. This
export is NOT data-free. It excludes input matrices, complete metadata, analysis
tables, model checkpoints and reports. It does not upload anything.
"""
from pathlib import Path
import argparse
import hashlib
from html.parser import HTMLParser
import json
import shutil
import urllib.parse
import zipfile


class Links(HTMLParser):
    def __init__(self):
        super().__init__();self.links=[]
    def handle_starttag(self, tag, attrs):
        if tag in {'a','img','script','link'}:
            self.links.extend(v for k,v in attrs if k in {'href','src'} and v)


def validate_site(root):
    paths=[p for p in root.rglob('*') if p.is_file()]
    if any(p.is_symlink() for p in root.rglob('*')):
        raise ValueError('Website must not contain symlinks')
    if not (root/'index.html').is_file():raise ValueError('Missing root index.html')
    for p in paths:
        if p.suffix.lower() not in {'.html','.txt','.json'} and p.name!='.nojekyll':
            raise ValueError('Unexpected website payload: '+str(p.relative_to(root)))
        if p.suffix=='.html':
            parser=Links();parser.feed(p.read_text())
            for link in parser.links:
                parts=urllib.parse.urlsplit(link)
                if parts.scheme or parts.netloc or not parts.path:continue
                target=(p.parent/urllib.parse.unquote(parts.path)).resolve()
                if root.resolve() not in target.parents and target!=root.resolve():
                    raise ValueError('Link leaves website: '+link)
                if not target.is_file():raise ValueError('Broken local link: '+link)
    size=sum(p.stat().st_size for p in paths)
    if size>=1_000_000_000:raise ValueError('Website reaches the 1 GB GitHub Pages site limit')
    return {'files':len(paths),'uncompressed_bytes':size,'input_matrices_included':False,
            'embedded_content':'Selected derived per-sample KO representation, rule evidence, coordinates, mapped environmental values, diagnostics and figure images are public when deployed.'}


def export(atlas, output):
    atlas=atlas.resolve();output=output.resolve()
    if output.exists():raise ValueError('Use a new output directory')
    manifest=json.loads((atlas/'results/viewer_manifest.json').read_text())
    if len(manifest['viewers'])!=8:raise ValueError('Expected all eight atlas panels')
    needed=['index.html']+list(manifest['viewers'].values())
    for name in needed:
        source=(atlas/name).resolve()
        if atlas not in source.parents or not source.is_file() or (atlas/name).is_symlink():
            raise ValueError('Invalid website source: '+name)
    output.mkdir(parents=True)
    for name in needed:
        dest=output/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(atlas/name,dest)
    shutil.copy2(atlas/'assets/vendor/leaflet-LICENSE.txt',output/'leaflet-LICENSE.txt')
    (output/'.nojekyll').write_text('')
    (output/'README.txt').write_text('Static mangrove functional atlas. Open index.html. This website embeds derived gene values, coordinates, environmental estimates and rendered images. Original matrices, full metadata, model outputs and checkpoints are not included. Detailed terrain requires internet. Source: https://github.com/fjbalvino/mangrove-sediment-biogeography\n')
    report=validate_site(output)
    report['profiles']=manifest['profiles']
    report['sha256']={p.relative_to(output).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.rglob('*')) if p.is_file()}
    (output/'site_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    archive=output.with_suffix('.zip')
    if archive.exists():raise ValueError('Archive already exists: '+str(archive))
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(output.rglob('*')):
            if p.is_file():z.write(p,p.relative_to(output))
    print(json.dumps({**report,'archive':str(archive),'archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest()},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--atlas',type=Path);parser.add_argument('--out',type=Path)
    parser.add_argument('--validate-site',type=Path)
    args=parser.parse_args()
    if args.validate_site:print(json.dumps(validate_site(args.validate_site.resolve()),indent=2))
    elif args.atlas and args.out:export(args.atlas,args.out)
    else:parser.error('Provide --atlas and --out, or --validate-site')
