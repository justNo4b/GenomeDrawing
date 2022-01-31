import argparse
from ssl import _create_unverified_context
from PIL import Image, ImageDraw, ImageFont
from collections import namedtuple
from argparse import ArgumentParser
from Bio import SeqIO

#define draw parameters
xSize       = 600
xName       = "Hospitalitermes Polycipiviridae-like virus"

padding_top         = 25
padding_left        = 45
padding_ORF         = 25
padding_ORF_level   = 15
padding_text        = 15
padding_text_x      = 10
padding_cap_y       = 37

ORF_width     = 10
line_width    = 5
isOnlyOne     = False

#define containers
workload     = []
segmentSizes = []
ORF_array    = []

##define ORF structure
ORF_struct = namedtuple("ORF_struct", "name start length drawlevel color")

#parse some arguments
parser = argparse.ArgumentParser()
parser.add_argument("-seg", "--file_input", nargs="*", type=str)
workload = parser.parse_args().file_input

#calculate some important stuff
seg_count   = len(workload)
thinkness   = padding_ORF / 2 + ORF_width + padding_ORF + padding_cap_y
#y-size с запасом, чтобы всё точно влезло
ySize       = int (xSize - padding_top * 4)


#####
# Work with seq
fileCount = 0
for file in workload:
    fileCount = fileCount + 1
    for sequence in SeqIO.parse(file, "gb"):
        segmentSizes.append(len(sequence.seq))
        allORFs  = []
        lastend = 0
        for j in range(0, len(sequence.features)):
            if(sequence.features[j].type == "CDS"):
                orfLevel = 0
                orfName  = ""
                orfColor = "Purple"
                orfStart = int(sequence.features[j].location.start)
                orfEnd   = int(sequence.features[j].location.end)
                if (lastend >= orfStart):
                    orfLevel = 1
                lastend  = orfEnd
                try:
                    orfLabel = sequence.features[j].qualifiers["label"][0]
                except:
                    orfLabel = ""
                try:
                    orfProd  = sequence.features[j].qualifiers["product"][0]
                except:
                    orfProd = ""

                if (orfLabel != ""):
                    orfName = orfLabel
                elif(orfProd != ""):
                    orfName = orfProd
                else:
                    orfName = "Unknown"
                if (orfProd == "RdRp" or orfProd == "polymerase"):
                    orfColor = "Green"
                newOrf = ORF_struct(orfName, orfStart , orfEnd - orfStart, orfLevel, orfColor)
                allORFs.append(newOrf)
        ORF_array.append(allORFs)

if (fileCount == 1):
    isOnlyOne = True

# найти масштаб, используя максимальную длину сгмента и xSize
# сперва найти самый длинный сегмент
max_length = 0
for i in range(0, len(segmentSizes)):
    if (segmentSizes[i] > max_length):
        max_length = segmentSizes[i]

# далее отмасштабировать, предполагая что мы хотим оставить по
# крайней мере padding справа и слева

scale = (xSize - padding_left * 2 ) / (max_length)
print(thinkness)

# теперь отрисовка
image = Image.new("RGBA", (xSize,ySize), (255,255,255,255))

# чтоб рисовать надо создать отдельный объект ImageDraw
drawObject = ImageDraw.Draw(image)
# создать шрифт для прорисовки
number_font    = ImageFont.truetype("OpenSans-Regular.ttf", 13)
header_font    = ImageFont.truetype("OpenSans-Regular.ttf", 16)
title_font    = ImageFont.truetype("OpenSans-Regular.ttf", 24)

nextSegY = padding_top
maxLevel = 0

# отрисовка циркулярных вирусов
padding_all = padding_top * 4
center_x    = (xSize) / 2
center_y    = (ySize) / 2
circle_xy   = (padding_all, padding_all, xSize - padding_all, ySize - padding_all)
circle_xy_a = []
circle_xy_a.append((padding_all + padding_ORF, padding_all + padding_ORF,  xSize - padding_all - padding_ORF, ySize - padding_all - padding_ORF))
circle_xy_a.append((padding_all + 2 * padding_ORF, padding_all + 2 * padding_ORF, xSize - padding_all - 2 * padding_ORF, ySize - padding_all - 2 * padding_ORF))

drawObject = ImageDraw.Draw(image)
drawObject.text((padding_left, padding_top), text=xName, fill="Black", font=title_font, anchor="la")
drawObject.ellipse(circle_xy, 'white', outline='Black', width=line_width)
# для циркулярных длина в центре круга
drawObject.text((center_x - padding_text_x, center_y - padding_text),  text = str(segmentSizes[i]), fill="Black", font=number_font, anchor="la")
#отрисовка ORF
c = i + 1
y = nextSegY + thinkness
x2 = padding_left + segmentSizes[i] * scale
for orf in ORF_array[i]:
    d_start = (orf.start * 360) / segmentSizes[i]
    d_end   = ((orf.start + orf.length) * 360) / segmentSizes[i]
    # отрисовка наших геномных блоков
    drawObject.arc(circle_xy_a[orf.drawlevel], d_start, d_end, orf.color, width=(ORF_width * 2))
    # отрисовка подписей
    #x_caption = int((x2_ORF + x1_ORF) / 2) - padding_text * 2
    #y_caption = int(y_ORF - ORF_width / 2)
    #drawObject.text((x_caption, y_caption),  orf.name, fill="White", font=number_font, anchor="la")
    nextSegY = nextSegY + padding_ORF_level * orf.drawlevel
nextSegY = nextSegY + thinkness


del drawObject
image.save("test2.png", "PNG")
