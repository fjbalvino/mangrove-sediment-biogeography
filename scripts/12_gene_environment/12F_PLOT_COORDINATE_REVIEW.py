#!/usr/bin/env python3
"""Plot externally reviewed coordinates and unresolved conflicts.

No geocoding is performed. Each map panel uses its actual ArcGIS EPSG:3857
export extent. Evidence and imagery remain external to the repository.
"""
import argparse
import json
import math
import io
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator
from PIL import Image

R = 6378137
def xy(lon, lat):
    return R*math.radians(lon), R*math.log(math.tan(math.pi/4+math.radians(lat)/2))
def lon(x):
    return math.degrees(x/R)
def lat(y):
    return math.degrees(2*math.atan(math.exp(y/R))-math.pi/2)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evidence", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    e = json.loads(a.evidence.read_text())
    base = a.evidence.parent
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11})
    fig, axs = plt.subplots(2,2,figsize=(16,15.5))
    fig.subplots_adjust(left=.065,right=.97,top=.88,bottom=.22,hspace=.47,wspace=.22)
    fig.suptitle(e["title"],x=.065,y=.97,ha="left",fontsize=23,fontweight="bold")
    fig.text(.065,.932,e["subtitle"],fontsize=13,color="#334155")
    for ax, panel in zip(axs.flat, e["panels"]):
        if panel["kind"] == "image":
            ax.imshow(Image.open(base/panel["image"]))
            ax.axis("off")
        else:
            export=json.loads((base/panel["export"]).read_text())
            extent=export["extent"]
            if extent["spatialReference"].get("latestWkid",extent["spatialReference"].get("wkid")) not in (3857,102100):
                raise ValueError("Map export must use EPSG:3857")
            im=Image.open(base/panel["image"])
            if im.size!=(export["width"],export["height"]):
                raise ValueError("Image size differs from export")
            ax.imshow(im,extent=[extent["xmin"],extent["xmax"],extent["ymin"],extent["ymax"]])
            x,y=xy(panel["lon"],panel["lat"])
            ax.scatter([x],[y],s=190,marker=panel.get("marker","o"),color=panel["color"],
                       edgecolors="white",linewidths=1.7,zorder=10)
            ax.annotate(panel["label"],(x,y),xytext=(12,20),textcoords="offset points",
                        color="white",fontsize=11,fontweight="bold",
                        bbox=dict(boxstyle="round,pad=.35",fc="#0f172a",ec="white",alpha=.9))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v,pos:f"{lon(v):.3f}° E"))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v,pos:f"{lat(v):.3f}° N"))
            ax.xaxis.set_major_locator(MaxNLocator(3));ax.yaxis.set_major_locator(MaxNLocator(3))
            ax.tick_params(labelsize=9)
            # Mercator distance adjusted to ground at panel latitude.
            dist=panel["scale_m"]/math.cos(math.radians(panel["lat"]))
            xx=extent["xmin"]+.06*(extent["xmax"]-extent["xmin"])
            yy=extent["ymin"]+.07*(extent["ymax"]-extent["ymin"])
            ax.plot([xx,xx+dist],[yy,yy],lw=4,color="white")
            ax.text(xx,yy+.025*(extent["ymax"]-extent["ymin"]),panel["scale_label"],
                    color="white",fontsize=10,bbox=dict(fc="#0f172a",alpha=.75,ec="none",pad=2))
        ax.set_title(panel["title"],loc="left",fontweight="bold",fontsize=13,pad=12)
        ax.text(0,-.12,panel["caption"],transform=ax.transAxes,ha="left",va="top",
                fontsize=10,color="#334155",linespacing=1.5)
    fig.text(.065,.04,e["footer"],fontsize=10,linespacing=1.5,color="#334155")
    buffer = io.BytesIO()
    fig.savefig(buffer,dpi=120,facecolor="white",format="png")
    buffer.seek(0)
    rendered = Image.open(buffer).convert("RGB")
    rendered.save(a.output,format="PNG",optimize=True)
    with Image.open(a.output) as check:
        check.load()
    plt.close(fig)

if __name__ == "__main__":
    main()

