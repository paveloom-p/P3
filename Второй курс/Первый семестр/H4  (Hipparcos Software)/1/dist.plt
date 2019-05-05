set terminal postscript eps 32                                                                                                  
set output 'distrib.eps'                                                                                                        
set nokey                                                                                                                       
set border 15 lw 3                                                                                                              
set size 1.5, 1                                                                                                                 
set xtics -8, 1, 16                                                                                                             
set label "     " at -13 , 3700  center rotate by 90                                                                            
set label "     " at -12 , 3800  center rotate by 90                                                                            
set label "     " at -11 , 3900  center rotate by 90                                                                            
set label "     " at -10 , 4000  center rotate by 90                                                                            
set label "     " at -9  , 4100  center rotate by 90                                                                            
set label "     " at -8  , 4200  center rotate by 90                                                                            
set label "2    " at -7  , 4302  center rotate by 90                                                                            
set label "8    " at -6  , 4408  center rotate by 90                                                                            
set label "30   " at -5  , 4530  center rotate by 90                                                                            
set label "100  " at -4  , 4700  center rotate by 90                                                                            
set label "310  " at -3  , 5010  center rotate by 90                                                                            
set label "1142 " at -2  , 5942  center rotate by 90                                                                            
set label "4233 " at -1  , 9133  center rotate by 90                                                                            
set label "11070" at  0  , 16070 center rotate by 90                                                                            
set label "19318" at  1  , 24218 center rotate by 90                                                                            
set label "17086" at  2  , 21886 center rotate by 90                                                                            
set label "15078" at  3  , 19778 center rotate by 90                                                                            
set label "10939" at  4  , 15539 center rotate by 90                                                                            
set label "5179 " at  5  , 9679  center rotate by 90                                                                            
set label "2357 " at  6  , 6757  center rotate by 90                                                                            
set label "1483 " at  7  , 5783  center rotate by 90                                                                            
set label "937  " at  8  , 5137  center rotate by 90                                                                            
set label "463  " at  9  , 4563  center rotate by 90                                                                            
set label "265  " at  10 , 4265  center rotate by 90                                                                            
set label "115  " at  11 , 4015  center rotate by 90                                                                            
set label "48   " at  12 , 3848  center rotate by 90                                                                            
set label "15   " at  13 , 3715  center rotate by 90                                                                            
set label "7    " at  14 , 3607  center rotate by 90                                                                            
set label "1    " at  15 , 3501  center rotate by 90                                                                            
plot [-8:16] [0:30000] 'result' u 1:2 w boxes lt -1 lw 10                                                                       
unset label                                                                                                                     
reset                                                                                                                           
