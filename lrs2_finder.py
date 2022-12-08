#simple code to make finder chart from PanSTARRS imaging
#  and overlay LRS2 IFU footprint on given coordinates with correct PA: lrs2_overlay.png
#   also outputs the PS2 r-filter FITS file cut-out on this position: ps2image.fits
#  
#   (python 3)
#   with much code copied from this workbook:
#     https://ps1images.stsci.edu/ps1_dr2_api.html
#
#  S. Janowiecki, Oct 2022



import sys
import re
import numpy as np
import matplotlib.pyplot as plt

import traceback
import argparse

import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord, FK5
from astropy.io.fits import writeto

 

def getimages(ra,dec,filters="grizy"):
    
    """Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = f"{service}?ra={ra}&dec={dec}&filters={filters}"
    table = Table.read(url, format='ascii')
    return table


def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
    
    """Get URL for images in the table
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """
    
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages(ra,dec,filters=filters)
    url = (f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           f"ra={ra}&dec={dec}&size={size}&format={format}")
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[np.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url





def getargs():
    myDescription = """
    Make LRS2 finder images from PS1 images (B and R have same angular orientation on sky)
    
    Examples:
        python lrs2_finder.py  -ra 00:04:35.1 -decneg Y -dec 00:15:52.3 -az E   (note negative dec!)
        python lrs2_finder.py  -ra 12:13:14.3 -decneg N -dec 15:16:10 -az E   
        python lrs2_finder.py  -ra 181.05503 -decneg N -dec 55.44353 -az W
         (either use decimal degrees for both, or sexagesimal for both - do not mix!)
    
    Image scaling parameter: percentile
        If image looks washed out and too white, use a larger value of pct, like 99.999
        If image looks too dark and sources are too faint, use a smaller value like 90


    """
    parser = argparse.ArgumentParser(description=myDescription,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)


    requiredNamed = parser.add_argument_group('required arguments')
    optNamed = parser.add_argument_group('optional arguments')
    requiredNamed.add_argument("-ra", help="RA of target, can be decimal degrees or sexagesimal", dest="ra_in", default="", required=True, type=str)
    requiredNamed.add_argument("-decneg", help="Is Dec negative? Y/N", dest="decneg_in", default="N", required=True, type=str)
    requiredNamed.add_argument("-dec", help="Dec of target, without sign, must match RA formatting", dest="dec_in", default="", required=True, type=str)
    #requiredNamed.add_argument("-ifu", help="IFU: B or R", dest="ifu_in", default="", required=True)
    requiredNamed.add_argument("-az", help="Azimuth: E or W", dest="az_in", default="", required=True)
    optNamed.add_argument("-pct", help="percentile for image scaling contrast", dest="pct", default=99.2, required=False)
    args = parser.parse_args()
    return args, parser





if __name__ == "__main__":

    #here we go
    ps1filename = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    fitscut = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
    
    (args, pars) = getargs()
    #print(args)
    

    #now parse azimuth request:
    if (args.az_in == 'E' or args.az_in == 'e'):
        az = 'E'
    elif (args.az_in == 'W' or args.az_in == 'w'):
        az = 'W'
    else:
        print('ERROR: Your AZ request is not recognized. Please select B or R')
        sys.exit()
    

    #convert into SkyCoord:
    #any colons?
    if (':' in args.ra_in or ':' in args.dec_in):
        #parsing sexagesimal:
        try:
            if (args.decneg_in=='Y' or args.decneg_in=='y'):
                c = SkyCoord(args.ra_in,'-'+args.dec_in, frame=FK5, unit=(u.hourangle, u.deg))
            else:
                c = SkyCoord(args.ra_in,args.dec_in, frame=FK5, unit=(u.hourangle, u.deg))
        except ValueError:
            traceback.print_exc()
            print('')
            print('Your angles are out of range, or there is a problem with negative Dec formatting: please try again')
            print('')
    else:
        #parsing decimal degrees
        try:
            if (args.decneg_in=='Y' or args.decneg_in=='y'):
                c = SkyCoord(args.ra_in, '-'+args.dec_in, unit="deg", frame=FK5)
            else:
                c = SkyCoord(args.ra_in, args.dec_in, unit="deg", frame=FK5)
        except ValueError:
            traceback.print_exc()
            print('')
            print('ERROR: Your angles are out of range, possibly a problem with negative Dec formatting, please try again')
            print('')

    #coordinate!
    #print(c)
    
    #now check for HET sanity:
    if (c.dec>71.5*u.deg or c.dec<-10.5*u.deg):
        print(' Your declination is outside of the range available to the HET (-10.5 to 71.5 degrees).')
        sys.exit()
    
    #and check if requesting W track that this declination gives two tracks:
    if (az == 'W'):
        if (c.dec>65.5*u.deg or c.dec<-4.3*u.deg):
            print('ERROR: you have requested a West track for a target at extreme Dec which has a single (Wast) track only')
            sys.exit()
    


    #whew, finally. if it gets here, then I think it's acceptable to continue!
    #time for trig from:  https://luna.mpe.mpg.de/wikihetdex/index.php/Sky_projection_with_HET
    ln = 30.681436*u.deg #deg
    lo = 55*u.deg #deg
    an = 0*u.deg #deg
    del0 = c.dec 

    #use some shortcuts to work in degrees
    sind = lambda degrees: np.sin(np.deg2rad(degrees))
    cosd = lambda degrees: np.cos(np.deg2rad(degrees))    
    asind = lambda rad: np.rad2deg(np.arcsin(rad))
    acosd = lambda rad: np.rad2deg(np.arccos(rad))
 
    #calculate azimuth:
    if (c.dec>65.5*u.deg):
        azdeg = 0*u.deg
        #pring(azdeg)
    elif (c.dec<4.3*u.deg):
        azdeg = 180*u.deg
        #print(azdeg)
    else:
        azdeg = acosd( (sind(del0) - (sind(ln)*sind(lo))) / (cosd(lo)*cosd(ln)) )
        azdeg2 = (180*u.deg - azdeg) + 180*u.deg
        #print(azdeg, azdeg2)

    #now for position angle:
    if (azdeg == 0*u.deg):
        padeg = 0*u.deg
        #print(padeg)
    elif (azdeg == 180*u.deg):
        padeg = 180*u.deg
        #print(padeg)
    else:
        padeg = acosd( ( (cosd(azdeg)*cosd(ln)*sind(lo)) - (cosd(lo)*sind(ln))) / (cosd(del0)) )
        padeg2 = (180*u.deg - padeg) + 180*u.deg
        #print(padeg,padeg2)
    
    #get correct Az & PA for this track
    if (az=='W'):
        azdeg = azdeg2
        padeg = padeg2
    #otherwise, keep the way it is
    
    print('az: '+str(azdeg))
    print('PA: '+str(padeg))
    #so this PA points in the direction of the short axis of the 12"x6" IFU
    

    #ok, wow, the hard part is done! (i hope) now just to get the PS1 image to match this
    #https://ps1images.stsci.edu/ps1image.html
    #following that example ^^
    from astropy.io import fits
    from astropy.visualization import PercentileInterval, AsinhStretch
    
    print('downloading FITS file from PanSTARRS...')
    print('')

    size=480 #pixels! these are 0.25" each, so this should give 120"x120"
    fitsurl = geturl(c.ra.deg, c.dec.deg, size=size, filters="r", format="fits")
    fh = fits.open(fitsurl[0])
    fim = fh[0].data
    # replace NaN values with zero for display
    fim[np.isnan(fim)] = 0.0
    #also inject fake value for image scaling
    fim[0,0]=9000000
    #image is flipped in Y direction!! ughhhh
    fim2 = np.flipud(fim)
    # set contrast to something reasonable
    transform = AsinhStretch() + PercentileInterval(float(args.pct))
    bfim = transform(fim2)    
    from scipy.ndimage import rotate
    bfim2 = rotate(bfim, 180*u.deg - padeg, axes=(1,0))

    #wcs = WCS(fh[0].header)
    plt.rcParams.update({'font.size':12})
    ax = plt.subplot()#projection=wcs)
    ax.imshow(bfim2,cmap="gray")



    ##label IFU position here
    xdim = np.shape(bfim2)[0]
    #print(xdim)
    ax.text(xdim/2,xdim/2,'         ', size=7,  va='center', ha='center',
            bbox=dict(boxstyle="square,pad=0.3", fc="none", ec="g", lw=1))
    #TESTED carefully to be sure it matches!
    #  these sizes and paddings are CRITICAL to getting this right. do NOT change!
    
    
    #helpful labels:
    ax.text(0,10,'RA: {:9.7f}deg'.format(c.ra.deg),size=8,color='w',va='center',ha='left')
    ax.text(0,30,'Dec: {:9.7f}deg'.format(c.dec.deg),size=8,color='w',va='center',ha='left')
    ax.text(0,50,'Track: '+str(az),size=8,color='w',va='center',ha='left')
    ax.text(0,70,'Az: {:6.3f}deg'.format(azdeg/u.deg),size=8,color='w',va='center',ha='left')
    ax.text(0,90,'PA: {:6.3f}deg'.format(padeg/u.deg),size=8,color='w',va='center',ha='left')
    

    ax.set_xticks([])
    ax.set_yticks([])



    plt.savefig('lrs2_overlay.png',bbox_inches='tight',dpi=200)
    
    #also save FITS file
    writeto('ps2image.fits',fh[0].data,header=fh[0].header,overwrite=True)


    
