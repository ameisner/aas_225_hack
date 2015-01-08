import os
import sqlite3

import lsst.afw.math as afwMath
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
from lsst.obs.sdss.processCcdSdss import ProcessCcdSdssTask
import tempfile
import numpy, math #cmath
import argparse
import Image

def scaleLinear(imageList, zmin=None, zmax=None,):
    '''Scale list of images using a linear ramp'''
    scaledImages = []

    for image in imageList:
        if (zmin==None):
            zmin = image.min()
        if (zmax==None):
            zmax = image.max()
        image = numpy.where(image > zmin, image, zmin)
        image = numpy.where(image < zmax, image, zmax)
        scaledImage = (image - zmin) * (255.0 / (zmax - zmin)) + 0.5

        # convert to 8 bit unsigned int ("b" in numpy)
        scaledImages.append(scaledImage.astype('uint8'))

    return scaledImages

def scaleAsinh(imageList, zmin=None, zmax=None, nonLinear = 3.0, scales=[4.9,5.7,7.8]):
    '''Scale list of images using asinh'''

    for image,scale in zip(imageList,scales):
        image -= image.min()
        image *= scale
        image = numpy.where(image < 20., image, 20.)

    # normalize images and apply asinh
    total = numpy.array(imageList).sum(axis=0)
    if (nonLinear == 0):
        value = total
    else:
        value = numpy.arcsinh(total*nonLinear)/nonLinear

    maxval = 1.
    for image in imageList:
        image *= value/total
        maxval = max(maxval, image.max())

    scaledImages = []
    for image in imageList:
        image /= maxval
        scaledImage = (image * 255.0 ) + 0.5
        # convert to 8 bit unsigned int ("b" in numpy)
        scaledImages.append(scaledImage.astype('uint8'))

    return scaledImages

def createRGB(imageList):
    '''Create RGB image from scale numpy arrays assuming order R,G,B'''
    mode='RGB'

    im = []
    for image in imageList:
        im.append(Image.fromarray(image, mode='L'))
    image = Image.merge('RGB', tuple(im))

    return image

def readFile(filename):
    '''Read exposure'''
    print filename
    exp = afwImage.ExposureF(filename)
    return exp

def makeCoord(ra, dec, coosys=afwGeom.degrees):
    '''Make coordinate'''
    return afwCoord.IcrsCoord(afwGeom.Angle(ra, coosys), afwGeom.Angle(dec, coosys))

def makeBbox(idx, dx, dy):
    '''Make bounding box'''
    xo = int(math.floor(idx.getX() - dx/2.))
    yo = int(math.floor(idx.getY() - dy/2.))
    return afwGeom.Box2I(afwGeom.Point2I(xo,yo), afwGeom.Extent2I(dx,dy))

def warpExposure(destExp, inpExp):
    '''warp exposure to match destExp'''
    warper = afwMath.Warper("lanczos2")
    omi = warper.warpImage(destExp.getWcs(), inpExp.getMaskedImage(),
                           inpExp.getWcs(), destBBox=destExp.getBBox(afwImage.PARENT))
    return omi


def coaddImages(ra, dec, expidlist, size, destDir, inputDir):
    '''Given an ra and dec coadd images given in list'''

    ra_dec = makeCoord(ra, dec)
    dx = dy = size
    offPt = afwGeom.Point2I(dx//2, dy//2)


    dataId = expidlist[0]
    #define first image as the template
    fullResult = ProcessCcdSdssTask.parseAndRun(
        args = [inputDir, "--output", destDir, "--id"] +
              ["%s=%s"%(key, val) for key, val in dataId.iteritems()],
	      doReturnResults = True,
    )
    butler = fullResult.parseCmd.butler
    bexp = butler.get('calexp', dataId)
    bwcs = bexp.getWcs()
    bbox = makeBbox(bwcs.skyToPixel(ra_dec)-offPt, dx, dy)

    try:
        bexp = bexp.Factory(bexp, bbox, True)
    except:
        raise ValueError("Bounding box does not overlap the initial image")

    #create coadded image
    bmiArr = bexp.getMaskedImage().getImage().getArray()
    bwcs = bexp.getWcs()
    nadded = 0
    for dataId in expidlist[1:]:
        fullResult = ProcessCcdSdssTask.parseAndRun(
           args = [inputDir, "--output", destDir, "--id"] +
	          ["%s=%s"%(key, val) for key, val in dataId.iteritems()],
	          doReturnResults = True,
        )
	butler = fullResult.parseCmd.butler
        exp = butler.get('calexp', dataId)
        twcs = exp.getWcs()
        bbox = makeBbox(twcs.skyToPixel(ra_dec), dx, dy)
        try:
            texp = exp.Factory(exp, bbox, True)
        except:
            print "Failed to create temp exposure from %s, probably because the bounding box is not in the image"%(file)
            continue
        nadded += 1
        wmi = warpExposure(bexp, texp)
        bmiArr += numpy.nan_to_num(wmi.getImage().getArray())

    print "Number of images added",nadded, bmiArr.max(), bmiArr.min()
    return bmiArr

def getExps(dataIds, nIm, inputDir, outPath):
    """inputDir: path to directory containing sql registry file
    """
    exps = {'g':[], 'r':[], 'i':[]}
    for dataId in dataIds:
        if nIm is None or len(exps[dataId['filter']]) < nIm:
            exps[dataId['filter']].append(dataId)
    return exps

def getDataIdsFromRaDec(ra, dec, sqlfile):
    conn = sqlite3.connect(sqlfile)
    result = conn.execute("select run, field, filter, camcol from fields"+\
                          " where ? between raMin and raMax and"+\
                          " ? between decMin and decMax")
    raise NotImplementedError("Not finished yet")


def main():
    '''Input a position on the sky and set of files (g,r,i) and return a color image'''
    # Setup command line options
    parser = argparse.ArgumentParser(description='Coadd R,G,B images and return color postage stamp')
    parser.add_argument('ra', type = float, help="RA")
    parser.add_argument('dec', type = float, help="dec")
    parser.add_argument('--outputFile', dest="outputFile", help='Output PNG file name', default='rgb.png')
    parser.add_argument('--destDir', dest="destDir", help='Destination directory.  Current directory by default.', default='.')
    parser.add_argument('--size', dest='size', type=int, default=300, help='Size of postage stamp')
    parser.add_argument('--nIm', dest='nIm', type=int, help='Number of images in each band to coadd.  Default is all of them')
    parser.add_argument('--sqlfile', dest='sqlfile', help='specify sqllite database file')
    parser.add_argument('--sqlregistrydir', dest='registryDir', help='path to directory containing registry file')
    args = parser.parse_args()

    outPath = tempfile.mkdtemp()
    sqlfile = args.sqlfile
    if not sqlfile:
        sqlfile = os.path.join(os.environ["HOME"], "AAS_Hack_Day/fields.sqlite")
    if not os.path.exists(sqlfile):
        raise RuntimeError("Could not find sql file: %s"%sqlfile)
    sqlregistrydir = args.sqlregistrydir
    if not sqlregistrydir:
        sqlregistrydir = os.path.join(os.environ["HOME"], "data/runs")
    if not os.path.exists(sqlregistrydir):
        raise RuntimeError("Could not find sql registry directory: %s"%sqlregistrydir)
    dataIds = getDataIdsFromRaDec(args.ra, args.dec, args.sqlfile)
    exps = getExps(dataIds, args.nIm, sqlregistrydir, outPath)
    # generate coadded images
    images = []
    for key, val in exps.iteritems:
        images.append(coaddImages(args.ra, args.dec, val, args.size), args.destDir)

    # scale image
#   scaledImages = scaleLinear(images, zmax=100.)
    scaledImages = scaleAsinh(images, zmax=100.)

    # create RGB
    image = createRGB(scaledImages)
    image.save(args.outputFile)

if __name__ == "__main__":
    main()
