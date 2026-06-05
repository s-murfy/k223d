#
#  This contains functions for
#          1) calculating an earthquake's size 
#          2) scaling relationships 
#
import numpy as np

def mag2area(Mw, author):
    # Provide estimate of rupture area based on magnitude and choosen empirical scaling law
    # input is moment magnitude (Mw) and author that scaling relationship is based on,
    # output is area in m**2
    
    kilo = 1000
    if author == "Wells&Coppersmith":
        area = 10.**(-3.99 + 0.98*Mw)  # this is for reverse faults 
    elif author == "Murotani":
        Mo = 10**(1.5*Mw+9.05)# convert from magnitude to moment 
        area = 1.34*10**(-10)*Mo**(2/3)
    elif author == "Skarlatoudis":  
        Mo = 10**(1.5*Mw+9.05)# convert from magnitude to moment 
        area =1.72*10**(-9)*Mo**0.62
    elif author ==  "Allen":
        if Mw <= 8.63:
            area = 10.**(-5.62+1.22*Mw)
        else:
            area = 10.**(2.23 +0.31*Mw)
    elif author ==  "Goda":
            area = 10**(-3.7135+0.9777*Mw)  #note: assuming no division between earthquake types 
    elif author ==  "AllenHayes_linear":
            area = 10**(-3.63+0.96*Mw)
    elif author ==  "AllenHayes_bilinear":
        if Mw<= 8.63:
            area = 10.**(-5.62+1.22*Mw)
        else:
            area = 10.**(2.23 +0.31*Mw)
    elif author ==  "Strasser":
        area = 10.**( -3.476 + 0.952*Mw)   #  scaling for interface events 
    else :
        print('author unrecognised, Strasser et al. 2010 assumed as default')
        area = 10.**( -3.476 + 0.952*Mw)         # default case is to use Strasser et al. 2010
        
    area = area*kilo**2    
    return area

def area2mag(area_m, author):
    # Provide estimate of rupture area based on magnitude and choosen empirical scaling law
    # input is area in m**2 and author that scaling relationship is based on,
    # output is moment magnitude (Mw
    
    kilo = 1000.
    area = area_m/kilo/kilo
    if author == "Wells&Coppersmith":
        mw = -1.61+0.41*np.log10(area)  # this is for reverse faults 
    elif author == "Murotani":
        Mo = (area/1.34*10**(10))**(3/2)
        Mw = (np.log10(Mo)-9.05)/1.5; # convert from magnitude to moment 
    elif author == "Skarlatoudis":    
        Mo = (area/1.72*10**(9))**(1./0.62)
        Mw = (np.log10(Mo)-9.05)/1.5; # convert from magnitude to moment 
    elif author ==  "Allen":
        if Mw <= 8.63:
            Mw = (5.62+np.log10(area))/1.22 # this is the re-ordering of area2mag equation
        else:
            Mw = (-2.23+np.log10(area))/0.31 # this is the re-ordering of area2mag equation
    elif author ==  "Goda":
            Mw = (3.713+np.log10(area))/0.9  #note: assuming no division between earthquake types 
    elif author ==  "AllenHayes_linear":
            Mw = (3.63+np.log10(area))/.96 # this is the re-ordering of area2mag equation
    elif author ==  "AllenHayes_bilinear":
        if Mw<= 8.63:
            Mw = (5.62+np.log10(area))/1.22 # this is the re-ordering of area2mag equation
        else:
            Mw = (-2.23+np.log(area))/0.31 # this is the re-ordering of area2mag equation
    elif author ==  "Strasser":
        Mw = 4.441+0.846*np.log10(area) #  scaling for interface events 
    else :
        print('author unrecognised, Strasser et al. 2010 assumed as default')
        Mw = 4.441+0.846*np.log10(area) #  scaling for interface events 
            
    return Mw
