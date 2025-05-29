# This code plots the output of the main C++ program. It is written in python
# since that way I can leverage matplotlib and the other excellent python
# ploting tools. Shapely is able to read WKT as well.

from pathlib import Path
from shapely import wkt
import matplotlib.pyplot as plt 
import geopandas as gp
plt.rcParams.update({'font.size': 22})
queryAttributes = {"facecolor":"k", "edgecolor":"k", "alpha":0.75}
neighborAttributes = {"facecolor":"b", "edgecolor":"k", "alpha":0.25}
titles = {"t":"Ground Truth for LSH Query", "q":"Results of LSH Query"}
labels = {"t":"", "q":"A"}

def jaccard_distance(g1, g2):
    return 1 - (g1.intersection(g2).area / g1.union(g2).area)

def plot_overlapping(filename):
    """
    Plot both the query polygon and the reported nearest neighbors in an
    overlapping fashion. The query will be a different color than the nearest
    neighbors.
    """
    output_filename = str(filename.with_name(f"{filename.stem}o.png"))

    fig =  plt.figure(figsize=(10,7.5))
    ax = fig.add_axes([1, 1, 1, 1])
    with open(str(filename), "r") as f:
        lines = f.readlines()
        polygons = [wkt.loads(l) for l in lines]
    s = gp.GeoSeries(polygons)
    gp.GeoSeries.plot(s.iloc[0:1], ax=ax, **queryAttributes)
    gp.GeoSeries.plot(s.iloc[1:], ax=ax, **neighborAttributes)
    ax.axis('off')
    ax.set_title(f"Overlapping {titles[filename.stem[-1]]} {int(filename.stem[0:-1])}")
    plt.savefig(output_filename, bbox_inches="tight")
    plt.close(fig)
    print(f"Plotted {filename.stem}o.png: {len(s)} polygons")


def plot_adjacent(filename):
    """
    Plot both the query polygon and the reported nearest neighbors with each
    one recieving a different plot. The query will be larger and in the center,
    with its neighbors around it.
    """
    output_filename = str(filename.with_name(f"{filename.stem}a.png"))

    fig =  plt.figure(figsize=(10,7.5))
    fig.suptitle(f"Split {titles[filename.stem[-1]]} {int(filename.stem[0:-1])}")
    with open(str(filename), "r") as f:
        lines = f.readlines()
        polygons = [wkt.loads(l) for l in lines]
    s = gp.GeoSeries(polygons)
    qax = fig.add_subplot(2, 2, 1)
    qax.set_title("Query Shape")
    qax.axis("off")
    gp.GeoSeries.plot(s.iloc[0:1], ax=qax, **queryAttributes)
    for i in range(1, len(s)):
        ax = fig.add_subplot(2, 2, i + 1, sharex=qax, sharey=qax)
        ax.set_title(f"{labels[filename.stem[-1]]}NN {i} (JD: {jaccard_distance(s[0], s[i]):.4})")
        ax.axis("off")
        gp.GeoSeries.plot(s.iloc[i:i+1], ax=ax, **neighborAttributes)
    fig.subplots_adjust(wspace=0)
    plt.savefig(output_filename, bbox_inches="tight")
    plt.close(fig)
    print(f"Plotted {filename.stem}a.png: {len(s)} polygons")


if __name__ == "__main__":
    files = [f for f in Path(r"../data/out/wkt/").rglob("*") if f.suffix == ".txt"]

    for file in files:
        plot_overlapping(file)
        plot_adjacent(file)

