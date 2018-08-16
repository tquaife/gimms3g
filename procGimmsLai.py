#!/usr/bin/env python

from datetime import date
from glob import glob
import os

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import matplotlib
import numpy as np
import scipy.stats as sp


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
    
    #whole years only
    for year in xrange(1982,2012):
      for month in monthList:
        for period in ["a","b"]:
          f=self.directory+"/AVHRRBUVI01."+str(year)+month+period+".abl"
          f=f.replace("//","/")
          if f in filenames:
            #print f
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

    return self.getTimeSeries_xy(x,y)


  def getTimeSeries_xy(self, x, y):

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

  
  def getDataCube(self,lat1,lon1,lat2,lon2):

    x1,y1 =self.getImgCoords(lat1,lon1)
    x2,y2 =self.getImgCoords(lat2,lon2)
  
    xMin=min([x1,x2])
    xMax=max([x1,x2])
    yMin=min([y1,y2])
    yMax=max([y1,y2])

    data=np.ones([len(self.fileList),yMax-yMin,xMax-xMin],dtype=np.uint8)

  
    for (n, fname) in enumerate(self.fileList):
      allData = np.fromfile(fname, dtype=np.uint8).reshape(self.ncols, self.nrows).T
      data[n,:,:]=allData[yMin:yMax,xMin:xMax]

    data=data.astype(np.float)*self.scale        
    data[np.where(data>(100*self.scale))]=np.nan    
        
    return data



  def getTrendImage(self,lat1,lon1,lat2,lon2):

    dataCube=self.getDataCube(lat1,lon1,lat2,lon2)  
    (t,yy,xx)=np.shape(dataCube)
    data=np.ones([yy,xx])*np.nan
    pimg=np.ones([yy,xx])*np.nan

    #divide by 24.0 to get in units lai per year
    points=np.arange(len(dataCube[:,0,0]))/24.0

    for x in xrange(xx):
      for y in xrange(yy):
        t=dataCube[:,y,x]
        idx = np.isfinite(t)
        if idx.sum()==0:
          continue
        #p=np.polyfit(points[idx],t[idx],1)
        #data[y,x]=p[0]
        (slope,intcpt,r,pval,stderr)=sp.linregress(points[idx],t[idx])
        data[y,x]=slope
        pimg[y,x]=pval
 
  
    return data, pimg



  
  def getImage(self, year, month, period):
    """Get one whole Earth image at the given date
    period should be "a" or "b".
    """
  
    fname = self.directory+"/AVHRRBUVI01."+str(year)+month+period+".abl"
    data = np.fromfile(fname, dtype=np.uint8).reshape(self.ncols, self.nrows)
    data = np.where(data>100, np.nan, data)
    data *= self.scale

    return data.T
 
 
 
  def plotNanHist(self,lat1,lon1,lat2,lon2):

    halfDate=self.dateList[int(len(self.dateList)/2.0)]
    hist1={}
    hist2={}
    
    x1,y1 =self.getImgCoords(lat1,lon1)
    x2,y2 =self.getImgCoords(lat2,lon2)
    
    xMin=min([x1,x2])
    xMax=max([x1,x2])
    yMin=min([y1,y2])
    yMax=max([y1,y2])
    
    for x in xrange(xMin,xMax): 
      for y in xrange(yMin,yMax): 
          
        ts=self.getTimeSeries_xy(x,y)
        
        #skip any always nan points
        if np.isfinite(ts).sum()==0:
          continue
        
        idx=np.where(np.isnan(ts))
        nanDates=np.array(self.dateList)[idx]
    
        for date in nanDates:
          hist=hist1
          if date > halfDate:
            hist=hist2
          idx=date.month*2
          if date.day==15:
            idx+=1
          if idx in hist:
            hist[idx]+=1
          else:
            hist[idx]=1
        
    print hist1
    print hist2

    hist3={}
    hist4={}
    for key in hist1:
      hist3[(key)/2.0]=hist1[key]
    for key in hist2:
      hist4[(key+0.25)/2.0]=hist2[key]

    plt.bar(hist3.keys(), hist3.values(), 0.3, color='g',alpha=0.7,label="1st half")
    plt.bar(hist4.keys(), hist4.values(), 0.3, color='r',alpha=0.7,label="2nd half")
    #plt.xlim([0,12])
    plt.xlabel("Month")
    plt.ylabel("Count")
    plt.legend()
    plt.show()


    
  def plotTimeSeriesTrend(self,lat,lon):
  
    ts=self.getTimeSeries(lat,lon)
    x=np.arange(len(ts))
    
    idx = np.isfinite(x) & np.isfinite(ts)
    poly=np.polyfit(x[idx],ts[idx],1)
    trend=np.polyval(poly,x)
    
    print poly
    plt.plot(x,ts)
    plt.plot(x,trend)
    plt.show()
  
 
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
    beg=mdates.date2num(self.dateList[0])
    end=mdates.date2num(self.dateList[-1])
    im=ax.imshow(hovmoller, interpolation=None,aspect='auto',\
                 extent=(beg,end,-90,90))
    date_format = mdates.DateFormatter('%Y')
    ax.xaxis_date()
    minor_locator = AutoMinorLocator(4)
    #minor_locator = MultipleLocator(365)
    ax.xaxis.set_minor_locator(minor_locator)
    cbar=plt.colorbar(im)
    cbar.set_label('LAI')
    ax.xaxis.grid(True,color='black',which='both')
    plt.title("Hovmoller plot of GIMMS 3g LAI for longitude "+str(lon))
    plt.xlabel("Date")
    plt.ylabel("Latitude")
    plt.tight_layout()
    if outfile==None:
      plt.show()      
    else:
      plt.savefig(outfile)
 


    

def doHovs(lonList):
  g=gimmsLai("/media/tqu/data/gimms")
  for lon in lonList:
    if lon<0:
      outf="hovmoller_%0.1fE.png"%abs(lon)
    else:
      outf="hovmoller_%0.1fW.png"%abs(lon)
    print outf
    g.hovmoller(lon,outfile=outf)





def doTrendPlots(lat1,lon1,lat2,lon2):

  from matplotlib.colors import ListedColormap
  from copy import copy
  
  g=gimmsLai("/media/tqu/data/gimms")
  
  (img,pvl)=g.getTrendImage(lat1,lon1,lat2,lon2)  
  
  pimg=copy(pvl)
  pimg[np.where(pvl<100.)]=3  
  pimg[np.where(pvl<0.10)]=2  
  pimg[np.where(pvl<0.05)]=1  
  pimg[np.where(pvl<0.01)]=0  
  
  f, axarr = plt.subplots(1, 2)
  
  im1=axarr[0].imshow(img,interpolation='none',cmap='Spectral',clim=(-0.04,0.04))
  axarr[0].set_title('Trend (LAI/year)')

  cMap = ListedColormap(['white', 'green', 'blue','red'])  
  im2=axarr[1].imshow(pvl,interpolation='none',cmap=cMap)
  axarr[1].set_title('p-vlaue')
  
  pos_ax=axarr[0].get_position()
  cax = plt.axes([pos_ax.xmin, 0.2, pos_ax.width, 0.03])
  plt.colorbar(im1,ax=axarr[0],cax=cax,orientation='horizontal')
  
  pos_ax=axarr[1].get_position()
  cax = plt.axes([pos_ax.xmin, 0.2, pos_ax.width, 0.03])
  tickPos=[0.5*(3./4.), 1.5*(3./4.), 2.5*(3./4.), 3.5*(3./4.)]
  cbar=plt.colorbar(im2,ax=axarr[1],cax=cax,orientation='horizontal',ticks=tickPos)
  
  cbar.ax.set_xticklabels(['p<0.01', '0.01<=p<0.05', '0.05<=p<0.1', '0.1<=p'])
  
  axarr[0].axis('off')
  axarr[1].axis('off')
  plt.show()
  




if __name__=="__main__":

 
  #doHovs([-55,-60,-65,20,10])

  doTrendPlots( 10.,-82.,-28.0,-40.)
  
  #g.plotNanHist(5,-75,-15,-54)
  #g.plotNanHist(0,-65,-1.5,-60.5)

  #g=gimmsLai("/media/tqu/data/gimms")
  #g.hovmoller(-65.,outfile="hovmoller_65E.png")
  #g.hovmoller(-65.)
  

  
