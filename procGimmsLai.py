#!/usr/bin/env python

from datetime import date
from glob import glob
import os

import matplotlib.pyplot as plt
import matplotlib
import numpy as np


monthList=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
day={'a':1,'b':15}

class gimmsLai(object):

  def __init__(self,directory="./"):
    """A basic class for handling GIMMS LAI data
    """
  
    self.nrows=2160
    self.ncols=4320
    self.resln=1./12.
    self.scale=0.1
    self.directory=directory

    self.genFileList()
  
  def genFileList(self):
    """Generate a list of files sorted in chronological order
    and a list of datetime date objects
    """

    self.fileList=[]
    self.dateList=[]
    filenames=glob(self.directory+"/*abl")
    
    for year in xrange(1980,2013):
      for month in monthList:
        for period in ["a","b"]:
          f=self.directory+"/AVHRRBUVI01."+str(year)+month+period+".abl"
          f.replace("//","/")
          if f in filenames:
            self.fileList.append(f)
            m=monthList.index(month)+1
            d=day[period]
            self.dateList.append(date(year=year,month=m,day=d))
 
  def getImgCoords(self, lat, lon):
    """Get the image corrdinates of a given
    latitude/longitude
    """
    y=int(np.floor((90.-lat)/self.resln))
    x=int(np.floor((180.+lon)/self.resln))

    return x, y
  
  
  def getTimeSeries(self, lat, lon):
    """Extract a time series of LAI data at
    a given latitude/longitude
    """
  
    x,y=self.getImgCoords(lat, lon)
  
    loc=x*self.nrows+y
  
    data=[]
    for fname in self.fileList:
    
      #do a seek to speed things up:
      f = open(fname, "rb")  
      f.seek(loc, os.SEEK_SET) 
      d=np.fromfile(f,dtype=np.uint8,count=1)[0]
      if d>100:
        data.append(np.nan)
      else:
        data.append(d*self.scale)
      
    return np.array(data)

  
  def getImage(self, year, month, period):
    """Get one whole Earth image at the given date
    period should be "a" or "b".
    """
  
    fname = self.directory+"/AVHRRBUVI01."+str(year)+month+period+".abl"
    data = np.fromfile(fname, dtype=np.uint8).reshape(self.ncols, self.nrows)
    data = np.where(data>100, np.nan, data)
    data *= self.scale

    return data.T
 
 
  def hovmoller(self,lon,outfile=None):
    """Show or save a hovmoller plot
    """ 
    
    x,y=self.getImgCoords(0., lon)
    hovmoller=np.zeros([self.nrows,len(self.dateList)])
    for (n,fname) in enumerate(self.fileList):
      #print fname
      data = np.fromfile(fname, dtype=np.uint8).reshape(self.ncols, self.nrows).T
      hovmoller[:,n]=data[:,x]
      
    hovmoller = np.where(hovmoller>100, np.nan, hovmoller)
    hovmoller *= self.scale
    
    matplotlib.rcParams.update({'font.size': 16})      
    fig, ax = plt.subplots(figsize=(18, 6))
    im=ax.imshow(hovmoller, interpolation=None,aspect=1/15.,extent=(1981,2012,-90,90))
    cbar=plt.colorbar(im)
    cbar.set_label('LAI')
    plt.title("Hovmoller plot of GIMMS 3g LAI for longitude "+str(lon))
    plt.xlabel("Date")
    plt.ylabel("Latitude")
    plt.tight_layout()
    if outfile==None:
      plt.show()      
    else:
      plt.savefig(outfile)
 


def test():

  year=2000
  month=7
  period="a"

  g=gimmsLai("/media/tqu/data/gimms")
  img=g.getImage(year,monthList[month-1],period)


  lat=0
  for lon in xrange(-90,160,10):
    x,y=g.getImgCoords(lat,lon)
    ts=g.getTimeSeries(lat,lon)
    i=g.dateList.index(date(year=year,month=month,day=day[period]))
  
    print g.dateList[i],lat,lon,x,y,ts[i],img[y,x]
    
  
  
if __name__=="__main__":


  g=gimmsLai("/media/tqu/data/gimms")
  g.hovmoller(-65.,outfile="hovmoller_65E.png")
  
  #test()

  
