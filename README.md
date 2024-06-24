<h1 align="center">GenomeDrawing</h1>

Small tool for in-scale drawing of the annotated viral genomes.
It takes a *.gb file(s) as an input and draw genome and annotated ORFs.

Requirements:

- PIL
- argparse
- SeqIO

Basic usage:

```
python3 draw.py -i input.gb -o output.png
```


Awaliable options:

- Basic

`-i` - input file or several files. If several *.gb entries are presented, it draws them as segments of the same virus

`-o` - Output file

- Drawing types

`-batch`    - if added when several input files are presented, draw them all as separate viruses

`-revcompl` - add this when sequence is a reverse-complement (draw 5` on the right etc)

`-circular` - if added draws all sequences a as if it they are circular

`-drawstops`- Draw in-frame stop codons (currently not compatible with circular mode)


- Visuals

`-width`    - (int) Width of the image (700 pixels as default)

`-name`     - (string) If drawn in segment mode this specifies a virus name

`-dpi`      - (int)   DPI of the obtained image, 300 is default

`-scale`    - (float) Scale everything bigger for high-res pictures
