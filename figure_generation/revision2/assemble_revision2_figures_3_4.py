#!/usr/bin/env python3
"""Assemble revision2 Figures 3 and 4 as manuscript TIFFs."""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


BUGFIX = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/"
    "codon2NS_manuscript/MBE/revision2"
)
FIGWORK = BUGFIX / "figwork"
OUT_DIR = BUGFIX / "submit/figs"


def load_font(size: int):
    for candidate in [
        Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
        Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf"),
    ]:
        if candidate.exists():
            return ImageFont.truetype(str(candidate), size=size)
    return ImageFont.load_default()


def open_rgb(path: Path) -> Image.Image:
    if not path.exists():
        raise FileNotFoundError(f"missing panel image: {path}")
    return Image.open(path).convert("RGB")


def fit_to_width(img: Image.Image, width: int) -> Image.Image:
    height = round(img.height * (width / img.width))
    return img.resize((width, height), Image.Resampling.LANCZOS)


def fit_square(img: Image.Image, size: int) -> Image.Image:
    return img.resize((size, size), Image.Resampling.LANCZOS)


def label(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, font) -> None:
    x, y = xy
    draw.text((x + 10, y - 28), text, fill="black", font=font)


def assemble_figure3() -> Image.Image:
    panel_width = 2460
    outer = 90
    gutter = 110
    font = load_font(170)

    a = fit_square(open_rgb(FIGWORK / "Figure_3A_codon_fitness_vs_log_expression_slope.png"), panel_width)
    b = fit_to_width(open_rgb(FIGWORK / "Figure_3B_observed_predicted_frequency_change_by_ghat_onepanel.png"), panel_width)

    width = outer * 2 + panel_width
    height = outer * 2 + a.height + gutter + b.height
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    pos_a = (outer, outer)
    pos_b = (outer, outer + a.height + gutter)
    canvas.paste(a, pos_a)
    canvas.paste(b, pos_b)
    label(draw, pos_a, "A", font)
    label(draw, pos_b, "B", font)
    return canvas


def assemble_figure4() -> Image.Image:
    panel_size = 2400
    outer = 90
    gutter = 110
    font = load_font(170)
    panels = [
        ("A", FIGWORK / "Figure_4A_factor_1_loadings_by_g.png"),
        ("B", FIGWORK / "Figure_4B_codon_stability_by_g.png"),
        ("C", FIGWORK / "Figure_4C_rna_stem_fold_change_by_2Ns.png"),
    ]

    width = outer * 2 + panel_size
    height = outer * 2 + panel_size * len(panels) + gutter * (len(panels) - 1)
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    for index, (panel_label, path) in enumerate(panels):
        img = fit_square(open_rgb(path), panel_size)
        pos = (outer, outer + index * (panel_size + gutter))
        canvas.paste(img, pos)
        label(draw, pos, panel_label, font)
    return canvas


def save(name: str, image: Image.Image) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tiff = OUT_DIR / f"Figure_{name}.tiff"
    preview = OUT_DIR / f"Figure_{name}_preview.png"
    image.save(tiff, format="TIFF", compression="tiff_lzw", dpi=(300, 300))
    image.save(preview, format="PNG", dpi=(300, 300))
    print(f"wrote {tiff}")
    print(f"wrote {preview}")


def main() -> None:
    save("3", assemble_figure3())
    save("4", assemble_figure4())


if __name__ == "__main__":
    main()
