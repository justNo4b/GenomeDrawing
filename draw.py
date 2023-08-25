import argparse
from PIL import Image, ImageDraw, ImageFont
from collections import namedtuple
from argparse import ArgumentParser
from Bio import SeqIO

#define draw parameters

padding_top         = 25
padding_left        = 45
padding_ORF         = 30
padding_ORF_level   = 50
padding_text        = 15
padding_text_x      = 10
padding_cap_y       = 37
ORF_width           = 13
line_width          = 5
isOnlyOne           = False
scale               = 1


def _genSeqName(seq):
    s = ""
    s += seq.annotations["organism"]
    return s

def _drawSegText(seqdata, batchmode, num):
    if (not batchmode):
        if (not isOnlyOne): drawObject.text((padding_left - padding_text, y - padding_cap_y), text=("Segment " + str(num)), fill="Black", font=header_font, anchor="la")
    # in batch mode treat each sequence as a separate virus
    # in that case pull name out of seq data
    else:
        drawObject.text((padding_left - padding_text, y - padding_cap_y), text=seqdata.name, fill="Black", font=header_font, anchor="la")
    return

def _defineStrandDirection(rcRNA):
    global leftGenomeText
    global rightGenomeText
    if (rcRNA == True):
        leftGenomeText  = "3`"
        rightGenomeText = "5`"
    else:
        leftGenomeText  = "5`"
        rightGenomeText = "3`"


def _drawComplementORF(drawObj, orf, currY, txtFont):
    y_ORF = currY + padding_ORF + padding_ORF_level * orf.drawlevel
    x1_ORF = padding_left + orf.start * scale
    x2_ORF = x1_ORF + orf.length * scale
    x_caption = int((x2_ORF + x1_ORF) / 2)
    y_caption = y_ORF
    x1_ORF_line = x1_ORF + ORF_width * 3
    x1_ORF_rect = x1_ORF + ORF_width * 2
    # отрисовка наших геномных блоков
    drawObj.rectangle([(x1_ORF_rect, y_ORF - ORF_width), (x2_ORF, y_ORF + ORF_width)], orf.color, outline="Black")
    drawObj.regular_polygon(bounding_circle=(x1_ORF_rect, y_ORF, ORF_width * 2), n_sides=3, rotation=90, fill=orf.color, outline="Black" )
    drawObj.line([(x1_ORF_line, y_ORF - ORF_width), (x1_ORF_line, y_ORF + ORF_width)], orf.color, width=1)
    # отрисовка подписей
    drawObj.text((x_caption, y_caption),  orf.name, fill="White", font=txtFont, anchor="mm")

    return

def _drawNormalORF(drawObj, orf, currY, txtFont):
    y_ORF = currY + padding_ORF + padding_ORF_level * orf.drawlevel
    x1_ORF = padding_left + orf.start * scale
    x2_ORF = x1_ORF + orf.length * scale
    x_caption = int((x2_ORF + x1_ORF) / 2)
    y_caption = y_ORF
    x2_ORF_line = x2_ORF - ORF_width * 3
    x2_ORF_rect = x2_ORF - ORF_width * 2
    # отрисовка наших геномных блоков
    drawObject.rectangle([(x1_ORF, y_ORF - ORF_width), (x2_ORF_rect, y_ORF + ORF_width)], orf.color, outline="Black")
    drawObj.regular_polygon(bounding_circle=(x2_ORF_rect, y_ORF, ORF_width * 2), n_sides=3, rotation=270, fill=orf.color, outline="Black" )
    drawObj.line([(x2_ORF_line, y_ORF - ORF_width), (x2_ORF_line, y_ORF + ORF_width)], orf.color, width=1)


    # отрисовка подписей
    drawObj.text((x_caption, y_caption),  orf.name, fill="White", font=txtFont, anchor="mm")

    return

#define containers
workload     = []
usedSequences = []
ORF_array    = []

##define structures
ORF_struct = namedtuple("ORF_struct", "name start length isCompl drawlevel color")
SEQ_data   = namedtuple("SEQ_data", "name length")

#parse some arguments
parser = argparse.ArgumentParser()
parser.add_argument("-seg", "--file_input", nargs="*", type=str)
parser.add_argument("-w", "--width", type=int, default=700)
parser.add_argument("-n", "--name", type=str, default="Virus name")
parser.add_argument("-o", "--output", type=str, default="default")
parser.add_argument("-rc", "--reverse_compl", action="store_true", help="Drawn image is reverse-complement")
parser.add_argument("-dpi", "--image_dpi", type=int, default=300, help="DPI of the obtained image, 300 is default")
parser.add_argument("-b", "--batch_mode", action="store_true", help="Treat each gb sequence as a separate virus")

workload    = parser.parse_args().file_input
xSize       = parser.parse_args().width
xName       = parser.parse_args().name
xOutput     = parser.parse_args().output
rcRNA       = parser.parse_args().reverse_compl
dpiUsed     = parser.parse_args().image_dpi
isBatch     = parser.parse_args().batch_mode


leftGenomeText  = ""
rightGenomeText = ""

#if outputfile dont specified and name is specified, set output name the same as name
if ((xOutput == "default") and (xName != "Virus name")):
    xOutput = xName


_defineStrandDirection(rcRNA)


#####
# Work with seq
seqCount = 0
for file in workload:
    for sequence in SeqIO.parse(file, "gb"):
        seqCount = seqCount + 1
        sLen = len(sequence.seq)
        sName = _genSeqName(sequence)
        newData = SEQ_data(sName, sLen)
        usedSequences.append(newData)
        allORFs  = []
        lastend = 0
        for j in range(0, len(sequence.features)):
            if(sequence.features[j].type == "CDS"):
                orfLevel    = 0
                orfName     = ""
                orfColor    = "Purple"
                orfStart    = int(sequence.features[j].location.start)
                orfEnd      = int(sequence.features[j].location.end)
                # Note -1 is complement, 1 is normal
                orfIsCompl = sequence.features[j].location.strand
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
                if (("RdRp" in orfProd) or ("polymerase" in orfProd)):
                    orfColor = "Green"
                newOrf = ORF_struct(orfName, orfStart , orfEnd - orfStart, orfIsCompl, orfLevel, orfColor)
                allORFs.append(newOrf)
        ORF_array.append(allORFs)

if (seqCount == 1):
    isOnlyOne = True


#calculate some important stuff
thinkness   = padding_ORF / 2 + ORF_width + padding_ORF + padding_cap_y
#y-size с запасом, чтобы всё точно влезло
ySize       = int (padding_top + thinkness * (seqCount + 1) * 1.5)
# найти масштаб, используя максимальную длину сгмента и xSize
# сперва найти самый длинный сегмент
max_length = 0
for i in range(0, len(usedSequences)):
    if (usedSequences[i].length > max_length):
        max_length = usedSequences[i].length

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
for i in range(0, len(usedSequences)):
    # calculate positions
    c = i + 1
    y = nextSegY + thinkness
    x2 = padding_left + usedSequences[i].length * scale
    #отрисовка общих вещей - геном, размер etc
    drawObject.line([(padding_left, y), (x2, y)], fill="Black", width=line_width)
    drawObject.text((padding_left - padding_text, y - padding_text), text=leftGenomeText, fill="Black", font=header_font, anchor="la")
    drawObject.text((int(x2 + padding_text / 2), y - padding_text), text=rightGenomeText, fill="Black", font=header_font, anchor="la")
    drawObject.text((x2, y + 5),  text = str(usedSequences[i].length), fill="Black", font=number_font, anchor="la")
    _drawSegText(usedSequences[i], isBatch, c)
    #отрисовка ORF
    for orf in ORF_array[i]:
        if (orf.isCompl == 1): _drawNormalORF(drawObject, orf, y, number_font)
        else: _drawComplementORF(drawObject, orf, y, number_font)
        nextSegY = nextSegY + padding_ORF_level * orf.drawlevel
    nextSegY = nextSegY + thinkness

# в конце не забыть удалить drawObject
del drawObject
image.save(xOutput, "PNG", dpi=(dpiUsed,dpiUsed))
