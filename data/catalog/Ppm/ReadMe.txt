         Catalogue of Positions and Proper Motions

                  S. Roeser and U. Bastian
         Astronomisches Rechen-Institut, Heidelberg
                         A. Kuzmin
          Sternberg Astronomical Institute, Moscow

Description:

  The PPM catalogue provides positions and proper motions for a total
  of 468861 stars, devided into four sections

    PPM North                   No.      1 - 181731
    PPM South                   No. 181732 - 378910
    PPM Bright Star Supplement  No. 400001 - 400321
    PPM Supplement              No. 700001 - 789676

  that have been combined into a single file PPM.dat. 

  PPM North gives J2000 positions and proper motions of 181,731
  stars north of -2.5 degrees declination. Its main purpose is to 
  provide a convenient, dense, and accurate net of astrometric 
  reference stars that represents the new IAU (1976) coordinate 
  system on the sky. The mean epoch is near 1931. The average mean 
  errors of the positions and proper motions are 0.27" and 0.43"/cen. 
  On the average six measured positions are available per star; 1064 
  stars do not have sufficient redundancies to resolve discrepancies.
  In addition to the positions and proper motions, the PPM (North)
  contains the BD number, the magnitude, the spectral type, the number
  of positions included, the mean error of each component of the
  position and proper motion, the weighted mean epoch in each
  coordinate, the numbers in the SAO, HD, and AGK3 catalogs, and 
  various standard notes.  In addition a number of individual notes 
  are given in the introduction. See PPMnorth.txt for more complete
  information.

  PPM South gives positions and proper motions of 197179
  stars south of about -2.5 degrees declination. Its main purpose
  is to provide a convenient, dense and accurate net of astrometric
  reference stars on the southern celestial hemisphere. This net is
  designed to represent as closely as possible the new IAU (1976)
  coordinate system on the sky, as defined by the FK5 star
  catalogue (Fricke et al., 1988).  In other words, it is a
   representation of this system at higher star densities and
  fainter magnitudes.  PPM South is the southern-hemisphere
  complement to the PPM Star Catalogue (Roeser and Bastian, 1991),
  which covers the northern hemisphere, plus a strip between the
  equator and about -2.5 degrees declination.  At the border line a
  continuous transition between the northern and southern part was
  tailored in such a way that no overlap, nor gaps, nor double
  entries occurred. See PPMsouth.txt for more complete information.

  A number of bright stars is missing from the PPM Star Catalogue, 
  both on the northern and on the southern hemisphere.
  The Bright Stars Supplement makes PPM complete down to V=7.5 mag.
  For this purpose it lists all missing stars brighter than
  V=7.6 mag that were found in published star lists. Their total
  number is 275. Only 2 of them are brighter than V=3.5. 
  The first version of the Bright Stars Supplement contained 321 stars,
  46 of which were erroneously included.  They are duplicates of stars
  already contained in the main parts of PPM. They are omitted here.
  Therefore the numbering sequence contains some gaps. See PPMbss.txt
  for more complete information.
  
  The 90,000 Stars Supplement to PPM is a subset of some 100,000
  CPC-2 stars from the Astrographic Catalogue that are not included
  in PPM South. See PPMsupp.txt for more complete information.


File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl     Records         Explanations
--------------------------------------------------------------------------------
ReadMe.txt        80           .         This file    
PPMnorth.txt      80        1145         Authors' description - PPM north 
PPMsouth.txt      80        1145         Authors' description - PPM south
PPMbss.txt        80        1145         Authors' description - PPM bright stars
PPMsupp.txt       80        1145         Authors' description - PPM supplement
PPM.dat          131      468861         Positions and Proper Motions Catalogue
--------------------------------------------------------------------------------


Byte-by-byte Description of file: PPM.dat
--------------------------------------------------------------------------------
   Bytes  Format  Units     Label    Explanations  
--------------------------------------------------------------------------------
   2-  7   I6     ---       PPM      Designation of the star
  10- 18   A9     ---       DM       BD (zones -02 to -22), else CD
  20- 23   F4.1   mag       mag      Magnitude
  25- 26   A2     ---       Sp       Spectral type
  28- 29   I2     h         RAh      Right Ascension J2000 (hours)
  31- 32   I2     min       RAm      Right Ascension J2000 (minutes)
  34- 39   F6.3   s         RAs      Right Ascension J2000 (seconds)
      42   A1     ---       DE-      Declination J2000 (sign)
  43- 44   I2     deg       DEd      Declination J2000 (degrees)
  46- 47   I2     arcmin    DEm      Declination J2000 (minutes)
  49- 53   F5.2   arcsec    DEs      Declination J2000 (seconds)
  56- 62   F7.4   s/yr      pmRA     Proper motion in RA, J2000
  64- 69   F6.3   arcsec/yr pmDE     Proper motion in DE, J2000
  71- 72   I2     ---       N        Number of individual published positions used
  74- 75   I2     10mas     e_RA     Mean error of RA
  77- 78   I2     10mas     e_DE     Mean error of DE
  80- 83   F4.1   mas/yr    e_pmRA   Mean error of pmRA
  85- 88   F4.1   mas/yr    e_pmDE   Mean error of pmDE
  90- 94   F5.2   yr        EpRA     Weighted mean epoch, RA and pmRA
  96-100   F5.2   yr        EpDE     Weighted mean epoch, DE and pmDE
 102-107   I6     ---       SAO      [1/258997]? SAO Designation
 109-114   I6     ---       HD       [1/350000]? Henry Draper Designation
 117-125   A9     ---       CPD      Cape Photographic Durchmusterung Designation
     127   A1     ---       Flag1    'P' - problem, 'C' - comment
     128   A1     ---       Flag2    'D' - double star
     129   A1     ---       Flag3    'H' - member of high-precision subset
     130   A1     ---       Flag4    'F' - member of FK5
     131   A1     ---       Flag5    'R' - remark, 'V' - V mag (CPC-2)
--------------------------------------------------------------------------------


