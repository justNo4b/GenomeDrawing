import argparse
from PIL import Image, ImageDraw, ImageFont
from collections import namedtuple
from argparse import ArgumentParser
from Bio import SeqIO

#define draw parameters

padding_top         = 25
padding_left        = 45
padding_ORF         = 25
padding_ORF_level   = 30
padding_text        = 15
padding_text_x      = 10
padding_cap_y       = 37

ORF_width     = 13
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
parser.add_argument("-w", "--width", type=int, default=700)
parser.add_argument("-n", "--name", type=str, default="Virus name")
parser.add_argument("-o", "--output", type=str, default="default")
workload = parser.parse_args().file_input
xSize = parser.parse_args().width
xName = parser.parse_args().name
xOutput = parser.parse_args().output
#if outputfile dont specified and name is specified, set output name the same as name
if ((xOutput == "default") and (xName != "Virus name")):
    xOutput = xName

#calculate some important stuff
seg_count   = len(workload)
thinkness   = padding_ORF / 2 + ORF_width + padding_ORF + padding_cap_y
#y-size с запасом, чтобы всё точно влезло
ySize       = int (padding_top + thinkness * (seg_count + 1) * 1.5)


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
title_font    = ImageFont.truetype("OpenSans-Regular.ttf", 18)

nextSegY = padding_top
maxLevel = 0

# Name on the top of it
drawObject.text((padding_left, padding_top), text=xName, fill="Black", font=title_font, anchor="la")
for i in range(0, len(segmentSizes)):
    # calculate positions
    c = i + 1
    y = nextSegY + thinkness
    x2 = padding_left + segmentSizes[i] * scale
    #отрисовка общих вещей - геном, размер etc
    drawObject.line([(padding_left, y), (x2, y)], fill="Black", width=line_width)
    drawObject.text((padding_left - padding_text, y - padding_text), text="5`", fill="Black", font=header_font, anchor="la")
    drawObject.text((int(x2 + padding_text / 2), y - padding_text), text="3`", fill="Black", font=header_font, anchor="la")
    if (not isOnlyOne):
        drawObject.text((padding_left - padding_text, y - padding_cap_y), text=("Segment " + str(c)), fill="Black", font=header_font, anchor="la")
    drawObject.text((x2, y + 5),  text = str(segmentSizes[i]), fill="Black", font=number_font, anchor="la")
    #отрисовка ORF
    for orf in ORF_array[i]:
        y_ORF = y + padding_ORF + padding_ORF_level * orf.drawlevel
        x1_ORF = padding_left + orf.start * scale
        x2_ORF = x1_ORF + orf.length * scale
        x_caption = int((x2_ORF + x1_ORF) / 2) - padding_text * 2
        y_caption = int(y_ORF - ORF_width / 2)
        # отрисовка наших геномных блоков
        drawObject.rectangle([(x1_ORF, y_ORF - ORF_width), (x2_ORF, y_ORF + ORF_width)], orf.color, outline="Black")
        # отрисовка подписей
        drawObject.text((x_caption, y_caption),  orf.name, fill="White", font=number_font, anchor="la")
        nextSegY = nextSegY + padding_ORF_level * orf.drawlevel
    nextSegY = nextSegY + thinkness

# в конце не забыть удалить drawObject
del drawObject
image.save(xOutput, "PNG")
