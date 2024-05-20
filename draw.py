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
scale_horizontal    = 1
leftGenomeText      = ""
rightGenomeText     = ""
num_font_size       = 13
head_font_size      = 16
title_font_size     = 18


def _genSeqName(seq):
    s = ""
    s += seq.annotations["organism"]
    return s


def _genLengthText(seqdata):
    s = ""
    if (seqdata.add5 != 0 or seqdata.add3 != 0):
        s += ">"
    s += str(seqdata.length - seqdata.add5 - seqdata.add3)
    return s

def _findStopReadThrough(sequence, start, end, list, add5, level):
    fragment = sequence.seq[int(start):int(end)].translate()
    size = 0
    tStr = fragment[:-1].split("*")
    for i in range(0, len(tStr) - 1):
        size += len(tStr[i]) * 3
        newStop = ORF_struct("STOP_CODON", start + add5 + size , 3, orfIsCompl, level, "", 0, 0, 0)
        list.append(newStop)
    return


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

def _drawComplementORF(drawObj, orf, currY, txtFont, drawArrowText, drawG3, drawG5):
    arrowColor = orf.color
    y_ORF = currY + padding_ORF + padding_ORF_level * orf.drawlevel
    x1_ORF = padding_left + orf.start * scale_horizontal
    x2_ORF = x1_ORF + orf.length * scale_horizontal
    x_caption = int((x2_ORF + x1_ORF) / 2)
    y_caption = y_ORF
    x1_ORF_line = x1_ORF + ORF_width * 3
    x1_ORF_rect = x1_ORF + ORF_width * 2
    if (drawG5 != 0): arrowColor = "Gray"
    # отрисовка наших геномных блоков
    if (drawArrowText):
        drawObj.rectangle([(x1_ORF_rect, y_ORF - ORF_width), (x2_ORF, y_ORF + ORF_width)], orf.color, outline="Black")
        drawObj.regular_polygon(bounding_circle=(x1_ORF_rect, y_ORF, ORF_width * 2), n_sides=3, rotation=90, fill=arrowColor, outline="Black" )
        drawObj.line([(x1_ORF_line, y_ORF - ORF_width), (x1_ORF_line, y_ORF + ORF_width)], orf.color, width=1)
        # отрисовка подписей
        drawObj.text((x_caption, y_caption),  orf.name, fill="White", font=txtFont, anchor="mm")
    else:
        drawObject.rectangle([(x1_ORF, y_ORF - ORF_width), (x2_ORF, y_ORF + ORF_width)], orf.color, outline="Black")

    return

def _drawNormalORF(drawObj, orf, currY, txtFont, drawArrowText, drawG3, drawG5):
    arrowColor = orf.color
    y_ORF = currY + padding_ORF + padding_ORF_level * orf.drawlevel
    x1_ORF = padding_left + orf.start * scale_horizontal
    x2_ORF = x1_ORF + orf.length * scale_horizontal
    x_caption = int((x2_ORF + x1_ORF) / 2)
    y_caption = y_ORF
    x2_ORF_line = x2_ORF - ORF_width * 3
    x2_ORF_rect = x2_ORF - ORF_width * 2
    if (drawG3 != 0): arrowColor = "Gray"
    # отрисовка наших геномных блоков
    # если регио слишком мелкий, опустить стрелку и подпись
    if(drawArrowText):
        drawObject.rectangle([(x1_ORF, y_ORF - ORF_width), (x2_ORF_rect, y_ORF + ORF_width)], orf.color, outline="Black")
        drawObj.regular_polygon(bounding_circle=(x2_ORF_rect, y_ORF, ORF_width * 2), n_sides=3, rotation=270, fill=arrowColor, outline="Black")
        drawObj.line([(x2_ORF_line, y_ORF - ORF_width), (x2_ORF_line, y_ORF + ORF_width)], orf.color, width=1)
        # отрисовка подписей
        drawObj.text((x_caption, y_caption),  orf.name, fill="White", font=txtFont, anchor="mm")
    else:
        drawObject.rectangle([(x1_ORF, y_ORF - ORF_width), (x2_ORF, y_ORF + ORF_width)], orf.color, outline="Black")

    return

def _drawORF(drawObj, orf, currY, txtFont, drawRC):
    x2_ORF = padding_left + (orf.start + orf.length) * scale_horizontal
    length = orf.length * scale_horizontal
    isTooSmall = length > ORF_width * 3
    print(str(isTooSmall) + " " + str(length) +" " + str(x2_ORF - ORF_width * 2))
    if (not drawRC): _drawNormalORF(drawObj, orf, currY, txtFont, isTooSmall, orf.drawG3, orf.drawG5)
    else: _drawComplementORF(drawObj, orf, currY, txtFont, isTooSmall, orf.drawG3, orf.drawG5)

    return

def _drawGap(drawObj, orf, currY):
    y_ORF = currY
    x1_ORF = padding_left + orf.start * scale_horizontal
    x2_ORF = x1_ORF + orf.length * scale_horizontal
    # отрисовка наших геномных блоков
    drawObj.rectangle([(x1_ORF, y_ORF - ORF_width), (x2_ORF, y_ORF + ORF_width)], orf.color, outline="Black")
    return

def _drawStop(drawObj, orf, currY, txtFont):
    y_ORF = currY
    x1_ORF = padding_left + orf.start * scale_horizontal
    x2_ORF = x1_ORF + orf.length * scale_horizontal
    # отрисовка наших геномных блоков
    drawObj.text((x1_ORF, y_ORF + ORF_width),  "*", fill="Black", font=txtFont, anchor="mm")
    return

#define containers
workload     = []
usedSequences = []
ORF_array    = []

##define structures
ORF_struct = namedtuple("ORF_struct", "name start length isCompl drawlevel color end drawG5 drawG3")
SEQ_data   = namedtuple("SEQ_data", "name length add5 add3")

#parse some arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--file_input", nargs="*", type=str, help="Output file")
parser.add_argument("-width", "--width", type=int, default=700, help="Width of the image (700 pixels as default)")
parser.add_argument("-name", "--name", type=str, default="Virus name", help="If drawn in segment mode this specifies a virus name")
parser.add_argument("-o", "--output", type=str, default="default.png", help="Output file")
parser.add_argument("-revcompl", "--reverse_compl", action="store_true", help="Sequence is a reverse-complement (draw 5` on the right)")
parser.add_argument("-dpi", "--image_dpi", type=int, default=300, help="DPI of the obtained image, 300 is default")
parser.add_argument("-batch", "--batch_mode", action="store_true", help="Treat each gb sequence as a separate virus")
parser.add_argument("-scale", "--scale_total", type=float, default=1.0, help="Scale everything bigger for high-res pictures")

workload    = parser.parse_args().file_input
xSize       = parser.parse_args().width
xName       = parser.parse_args().name
xOutput     = parser.parse_args().output
rcRNA       = parser.parse_args().reverse_compl
dpiUsed     = parser.parse_args().image_dpi
isBatch     = parser.parse_args().batch_mode
scale       = parser.parse_args().scale_total





#if outputfile dont specified and name is specified, set output name the same as name
if ((xOutput == "default") and (xName != "Virus name")):
    xOutput = xName


_defineStrandDirection(rcRNA)
padding_top         = int( scale * padding_top)
padding_left        = int( scale * padding_left)
padding_ORF         = int( scale * padding_ORF)
padding_ORF_level   = int( scale * padding_ORF_level)
padding_text        = int( scale * padding_text)
padding_text_x      = int( scale * padding_text_x)
padding_cap_y       = int( scale * padding_cap_y)
ORF_width           = int( scale * ORF_width)
line_width          = int( scale * line_width)
num_font_size       = int( scale * num_font_size)
head_font_size      = int( scale * head_font_size)
title_font_size     = int( scale * title_font_size)

#####
# Work with seq
seqCount = 0
for file in workload:
    for sequence in SeqIO.parse(file, "gb"):
        add5 = 0
        add3 = 0
        seqCount = seqCount + 1
        sLen = len(sequence.seq)
        sName = _genSeqName(sequence)
        allORFs  = []
        lastend = 0
        orfLevel    = 0
        for j in range(0, len(sequence.features)):
            if(sequence.features[j].type == "CDS" or sequence.features[j].type == "gap"):

                orfName     = ""
                orfColor    = "Purple"
                orfStart    = int(sequence.features[j].location.start)
                orfEnd      = int(sequence.features[j].location.end)
                # Note -1 is complement, 1 is normal
                orfIsCompl = sequence.features[j].location.strand
                if (lastend > orfStart):
                    orfLevel += 1
                else: orfLevel = 0
                lastend  = orfEnd

                if (sequence.features[j].type == "CDS"):
                    _findStopReadThrough(sequence, orfStart, orfEnd, allORFs, add5, orfLevel)

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
                if (("RdRp" in orfProd) or ("polymerase" in orfProd) or (("RdRp" in orfLabel))):
                    orfColor = "Green"

                ################################################################
                #   Treatment of the gaps
                #
                if (sequence.features[j].type == "gap"):
                    orfName  = "GAP"
                    orfColor = "Gray"
                    print(str(orfStart) + "    .. " + str(sLen - add5))
                    if (orfStart <= 1 and orfEnd == 1):
                        add5 = int(sequence.features[j].qualifiers["estimated_length"][0])
                        sLen += add5
                        orfStart = -add5
                        orfEnd = 0
                    elif (orfStart == (sLen - add5 - 1) and orfEnd == (sLen - add5)):
                        add3 = int(sequence.features[j].qualifiers["estimated_length"][0])
                        #####################
                        # 3-end
                        # Loop all existing features. If ORFs finish 3 nt to the end, extend them by add3 as well
                        #
                        for i in range(0, len(allORFs)):
                            orf = allORFs[i]
                            if (orf.end >= sLen - add5 - 2):
                                #gappedORF = ORF_struct("", orf.end, add3 + 1, 0, orf.drawlevel, "Gray", -1)
                                allORFs[i] = allORFs[i]._replace(length = (sLen + add3 + 1 - orf.start))
                                allORFs[i] = allORFs[i]._replace(drawG3 = add3)
                        orfStart += 1
                        sLen += add3
                        orfEnd = sLen - add5 + 1
                        #allORFs.insert(0, gappedORF)


                newOrf = ORF_struct(orfName, orfStart + add5 , (orfEnd + add5) - (orfStart + add5), orfIsCompl, orfLevel, orfColor, (orfEnd + add5), 0, 0)
                allORFs.append(newOrf)
        newData = SEQ_data(sName, sLen, add5, add3)
        usedSequences.append(newData)
        ORF_array.append(allORFs)

if (seqCount == 1):
    isOnlyOne = True


#calculate some important stuff
thinkness   = padding_ORF / 2 + ORF_width + padding_ORF + padding_cap_y
#y-size с запасом, чтобы всё точно влезло
ySize       = int (padding_top + thinkness * (seqCount + 1) * 1.5 * scale)
# найти масштаб, используя максимальную длину сгмента и xSize
# сперва найти самый длинный сегмент
max_length = 0
for i in range(0, len(usedSequences)):
    if (usedSequences[i].length > max_length):
        max_length = usedSequences[i].length

# далее отмасштабировать, предполагая что мы хотим оставить по
# крайней мере padding справа и слева

scale_horizontal = (xSize - padding_left * 2 ) / (max_length)

# теперь отрисовка
image = Image.new("RGBA", (xSize,ySize), (255,255,255,255))

# чтоб рисовать надо создать отдельный объект ImageDraw
drawObject = ImageDraw.Draw(image)
# создать шрифт для прорисовки
number_font    = ImageFont.truetype("OpenSans-Regular.ttf", num_font_size)
header_font    = ImageFont.truetype("OpenSans-Regular.ttf", head_font_size)
title_font    = ImageFont.truetype("OpenSans-Regular.ttf", title_font_size)

nextSegY = padding_top
maxLevel = 0

# Name on the top of it
drawObject.text((padding_left, padding_top), text=xName, fill="Black", font=title_font, anchor="la")
for i in range(0, len(usedSequences)):
    # calculate positions
    c = i + 1
    y = nextSegY + thinkness
    x2 = padding_left + usedSequences[i].length * scale_horizontal
    #отрисовка общих вещей - геном, размер etc
    drawObject.line([(padding_left, y), (x2, y)], fill="Black", width=line_width)
    drawObject.text((padding_left - padding_text, y - padding_text), text=leftGenomeText, fill="Black", font=header_font, anchor="la")
    drawObject.text((int(x2 + padding_text / 2), y - padding_text), text=rightGenomeText, fill="Black", font=header_font, anchor="la")
    drawObject.text((x2 + 5, y + 5),  text = _genLengthText(usedSequences[i]), fill="Black", font=number_font, anchor="la")
    _drawSegText(usedSequences[i], isBatch, c)
    #отрисовка ORF
    for orf in ORF_array[i]:
        if (orf.name == "GAP"):
            _drawGap(drawObject, orf, y)
        elif (orf.name == "STOP_CODON"):
            _drawStop(drawObject, orf, y, number_font)
        else:
            if (orf.isCompl == 1): _drawORF(drawObject, orf, y, number_font, False)
            else: _drawORF(drawObject, orf, y, number_font, True)
        nextSegY = nextSegY + padding_ORF_level * orf.drawlevel
    nextSegY = nextSegY + thinkness

# в конце не забыть удалить drawObject
del drawObject
image.save(xOutput, "PNG", dpi=(dpiUsed,dpiUsed))
