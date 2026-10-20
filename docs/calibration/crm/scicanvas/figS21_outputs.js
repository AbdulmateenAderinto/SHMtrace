() => {
  // SHMtrace's outputs at a glance: what a project run writes, and how to read the figures. The top band is the
  // run (inputs, SHMtrace, four kinds of output); below, the 26 figures grouped by the question each answers,
  // each group with one of its figures redrawn from the J_H4 example run and a short note on how to read it.
  // Numbers are written into DATA by outputs_data.py. Built with build.js into crm/figures/FigureS21_outputs;
  // the six question cards are panels A-F.
  // Cell Press rules as in crm_figure1.js (174 mm; Arial 6-8 pt; strokes 0.5-1.5 pt). Wild type blue,
  // AID-null red (also the controls' colour, as in every SHMtrace figure).
  const DATA = /*DATA*/{"qc":{"WT2":[100.0,68.28,60.8,56.02,10.0],"WT3":[100.0,67.76,46.04,43.28,10.0],"WT4":[100.0,62.41772707327775,24.023694602896008,21.105748135146996,10.969723562966212],"KO1":[100.0,64.42,32.72,29.98,10.0],"KO2":[100.0,56.66,27.26,23.36,10.0],"KO3":[100.0,52.04,30.54,24.88,10.0]},"bg":{"pos":[101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629],"rate":[0.00080766,0.00127006,0.00195965,0.00057571,0.00057544,0.00103603,0.00103794,0.00034566,0.00126714,0.00103699,0.00057544,0.00080562,0.00080562,0.00057544,0.00057544,0.00057571,0.00011514,0.00034534,0.00011509,0.00011535,0.00034503,0.00011506,0.00034622,0.00057504,0.00057557,0.00103532,0.00034582,0.00034566,0.0012686,0.00196146,0.00034542,0.00287389,0.00126451,0.00034487,0.00057478,0.00057491,0.00034495,0.00011501,0.00011498,0.00034495,0.00011496,0.00034487,0.00563413,0.00034495,0.00011498,0.00011498,0.00011498,0.00103508,0.00080488,0.00011498,0.00057478,0.00034622,0.00080562,0.00034542,0.00057531,0.00080654,0.00034495,0.00080543,0.00034503,0.00103484,0.00057491,0.00011496,0.00034495,0.00034487,0.00080469,0.00080506,0.00057504,0.0010346,0.00126509,0.00241518,0.00126509,0.00034126,0.00420885,0.00034141,0.00011375,0.0003411,0.00056863,0.00034149,0.00056876,0.00034134,0.00193512,0.00011375,0.00147946,0.00102354,0.00056863,0.00079645,0.0019338,0.00445155,0.00056902,0.001251,0.00034141,0.00079954,0.00056941,0.00079954,0.00011375,0.00147879,0.00079645,0.00102401,0.00057071,0.00284317,0.00970874,0.00238772,0.00056876,0.00079627,0.00420115,0.00063028,0.00652996,0.00034204,0.00466917,0.00034227,0.00125156,0.00079573,0.00056863,0.00329508,0.01384915,0.00034087,0.00466068,0.00263248,0.00056863,0.00443232,0.00170474,0.00238718,0.00102284,0.00124986,0.00147744,0.00034118,0.00147811,0.00011375,0.00147879,0.00034118,0.00080137,0.00125556,0.00056902,0.00170474,0.00102261,0.00034095,0.00238663,0.00034095,0.00125043,0.00056799,0.00056838,0.00034072,0.001251,0.00056786,0.000795,0.00034775,0.00011409,0.00443232,0.00216229,0.00261572,0.0005685,0.00102471,0.00193116,0.00034087,0.0014771,0.0017059,0.00125014,0.00261334,0.00306644,0.00034072,0.00124957,0.000795,0.00034118,0.00011362,0.00034079,0.00034079,0.00193072,0.00034079,0.00124901,0.00083403,0.00034064,0.00124901,0.00102261,0.00080562,0.0007979,0.00034056,0.00192984,0.00011352,0.00102191,0.00034064,0.00079554,0.00011357,0.0005676,0.00034056,0.00056747,0.00011349,0.00079482,0.00034048,0.00011347,0.00011344,0.00034033,0.00056721,0.00034033,0.00011344,0.00034041,0.00034041,0.00011347,0.00034048,0.00034048,0.00011365,0.00147643,0.00625071,0.00034087,0.00102261,0.00056863,0.00170358,0.00079464,0.00102215,0.00124872,0.00147879,0.00034048,0.00170203,0.00056721,0.00056863,0.00102076,0.00011342,0.00056709,0.00102099,0.00034041,0.00056721,0.00011378,0.00011342,0.00011342,0.00034025,0.00034025,0.00011365,0.00079428,0.00011342,0.00079392,0.00056709,0.00102076,0.00056747,0.0007941,0.00124759,0.00056786,0.0007941,0.00011347,0.00056812,0.00102215,0.00238447,0.00034048,0.01667612,0.00079428,0.00102471,0.00034025,0.00079374,0.00034017,0.00011339,0.00034165,0.00215885,0.00034087,0.00011339,0.00056747,0.00034025,0.00056696,0.00034017,0.0007941,0.00034017,0.00034079,0.00056721,0.00079374,0.00034017,0.00034017,0.00079392,0.00034017,0.00011339,0.00124787,0.00147442,0.00079428,0.00056721,0.00034033,0.00034033,0.00011344,0.00034025,0.00011339,0.00034041,0.00079392,0.00056696,0.00034033,0.00011339,0.00056709,0.00011342,0.00079356,0.00079392,0.00056696,0.00079428,0.00283543,0.00215885,0.00034025,0.00034002,0.00079338,0.00079356,0.00034056,0.00034056,0.0003401,0.00147376,0.00079338,0.0007941,0.00079374,0.00056709,0.00192766,0.00011339,0.0007941,0.00034017,0.00056889,0.00079338,0.00034002,0.00147442,0.00308184,0.00011344,0.00056683,0.00034376,0.00057097,0.00056683,0.00034002,0.00102006,0.00079374,0.00192809,0.00056799,0.00283479,0.00147409,0.0005676,0.00034017,0.00011339,0.00011334,0.00170087,0.00034033,0.0005667,0.00034002,0.00079338,0.0003401,0.00056683,0.00124702,0.00011334,0.00011339,0.00034056,0.00147409,0.00056683,0.00034041,0.0003401,0.00102052,0.00102029,0.00079392,0.00260801,0.00102076,0.0003429,0.00079464,0.01082621,0.00034017,0.00056747,0.00056709,0.00034033,0.00079374,0.00056683,0.00011337,0.00011334,0.00034002,0.00034025,0.0003401,0.00011337,0.00034002,0.0003401,0.00034002,0.0003401,0.0005667,0.00011334,0.00056683,0.00056683,0.00034025,0.00034033,0.0005667,0.0003401,0.00102029,0.00034002,0.00056747,0.00079392,0.00056683,0.00079392,0.00011339,0.00056696,0.0003401,0.00079392,0.00170087,0.00397231,0.00011339,0.00056683,0.00102191,0.0003401,0.00056683,0.00056683,0.00170358,0.00147543,0.0003401,0.00079356,0.00034025,0.0003401,0.00034103,0.00011342,0.00056683,0.00079464,0.00056683,0.00056683,0.00011339,0.00170396,0.00102076,0.00034002,0.00102099,0.0005667,0.0005685,0.00011349,0.00056683,0.00079374,0.00079392,0.00147376,0.00147409,0.00011344,0.00056683,0.00102029,0.00352073,0.00148487,0.00056709,0.00056721,0.00056747,0.00102029,0.0003401,0.00056734,0.00011337,0.00011339,0.00034017,0.0005667,0.00124759,0.00056709,0.00238176,0.00056825,0.00307622,0.00011334,0.00124702,0.00147342,0.00125156,0.00487916,0.00534638,0.00443636,0.00056825,0.00124702,0.00011344,0.00623795,0.00510146,0.00079392,0.00124731,0.00034025,0.0003401,0.00124731,0.00192809,0.00056747,0.00056696,0.00079338,0.00034033,0.00056696,0.0005667,0.0005667,0.00034002,0.00034017,0.00079374,0.00102029,0.00011337,0.00011342,0.00011337,0.00011337,0.00011342,0.00034033,0.00034017,0.00034056,0.00034025,0.00011342,0.00034025,0.00034072,0.00034064,0.00034048,0.0003411,0.00193028,0.00170435,0.00034056,0.00011349,0.00034048,0.00011368,0.00079464,0.0001137,0.00102168,0.00034056,0.00011349,0.00034056,0.00079482,0.0005676,0.00079446,0.00034048,0.00147845,0.00102308,0.00102191,0.00011357,0.00216278,0.00239918,0.01148641,0.00125642,0.0003401,0.00170087,0.00192897,0.00923287,0.00830584,0.00124731,0.00124674,0.0005667,0.0005667,0.00056683,0.0017001,0.00056696,0.00193072,0.00102145,0.00124674,0.00102099,0.00147543,0.00079374,0.00056734,0.00215885,0.00306505,0.0007941,0.00102099,0.00056721,0.00056747,0.00102029]},"dens":{"rows":[["WT2","Wt",83.748,0.378952],["WT3","Wt",133.205,0.387223],["WT4","Wt",27.826,0.381177],["KO1","Aicda-/-",0.0,0.379331],["KO2","Aicda-/-",0.76,0.379856],["KO3","Aicda-/-",0.0,0.380227]],"p":0.002672,"fold":141.496427},"spec":{"classes":["A>T","A>C","A>G","T>A","T>C","T>G","C>A","C>T","C>G","G>A","G>T","G>C"],"group":[9.561128526645769,9.404388714733543,18.495297805642632,8.463949843260188,14.576802507836991,5.799373040752351,3.134796238244514,7.6802507836990594,4.231974921630094,8.77742946708464,4.545454545454546,5.329153605015674],"animals":[[9.049773755656108,9.95475113122172,17.647058823529413,11.312217194570136,13.122171945701357,5.882352941176471,2.7149321266968327,6.787330316742081,3.6199095022624435,9.95475113122172,5.429864253393665,4.524886877828054],[10.755813953488373,9.883720930232558,19.476744186046513,7.267441860465116,14.244186046511627,4.941860465116279,3.7790697674418605,7.267441860465116,4.069767441860465,8.430232558139535,4.3604651162790695,5.523255813953488],[5.47945205479452,5.47945205479452,16.438356164383563,5.47945205479452,20.54794520547945,9.58904109589041,1.36986301369863,12.32876712328767,6.8493150684931505,6.8493150684931505,2.73972602739726,6.8493150684931505]],"gc_hot":[56.2,48.7,55.6]},"pos":{"wt":[[101,0.5333],[102,0.1333],[105,0.3333],[107,7.8],[108,1.8],[109,3.3333],[110,1.4667],[111,6.6667],[112,1.0],[114,0.6],[117,1.0],[121,0.3333],[122,0.8],[123,0.4667],[124,0.4],[125,1.6],[127,0.7333],[128,0.2],[129,0.6667],[130,0.2],[131,0.1333],[132,1.5333],[133,1.6],[135,1.6],[136,2.7333],[137,2.2],[138,5.8667],[139,4.0],[140,0.1333],[141,2.7333],[143,9.1333],[144,1.2667],[145,0.8667],[146,4.8667],[147,1.2667],[149,0.3333],[151,7.8],[152,1.8],[154,0.5333],[156,2.9333],[157,1.3333],[159,0.6667],[160,15.9333],[161,0.8667],[162,0.6],[163,0.2667],[164,1.6667],[165,7.0667],[166,11.6667],[167,3.2],[168,0.8667],[169,0.0667],[170,2.8],[171,1.0],[172,3.6],[173,1.9333],[175,1.0667],[176,1.0],[177,1.0],[178,1.8],[179,0.2667],[180,0.2],[181,1.1333],[182,2.6667],[183,3.0667],[184,1.7333],[185,0.0667],[186,0.1333],[187,0.0667],[188,0.8667],[191,1.7333],[192,2.0],[193,0.4667],[197,2.4667],[198,1.0],[199,1.5333],[202,2.3333],[206,0.7333],[207,0.5333],[208,2.0],[209,2.6667],[210,0.3333],[211,1.8],[212,6.4],[213,1.3333],[214,1.2],[216,1.0667],[217,3.8],[218,2.4667],[219,5.4667],[220,0.1333],[221,3.2],[222,3.8],[223,2.9333],[224,2.0],[225,1.6667],[228,0.2],[230,1.6],[231,1.6667],[235,0.8],[237,0.8667],[240,0.0667],[242,0.1333],[244,1.2],[245,0.4],[246,0.0667],[247,2.0],[248,0.2667],[252,2.3333],[253,1.5333],[255,0.0667],[257,0.1333],[258,3.9333],[259,1.4],[261,1.4],[263,1.1333],[264,0.4667],[266,0.5333],[267,0.1333],[269,1.5333],[270,0.5333],[271,0.9333],[273,1.7333],[274,1.3333],[275,4.3333],[276,0.0667],[277,1.1333],[278,0.8667],[279,7.4667],[280,0.2667],[281,4.2],[282,3.8],[283,0.2667],[284,1.2],[285,6.8],[286,1.1333],[287,3.1333],[289,0.7333],[290,1.1333],[291,0.1333],[293,1.0],[294,3.2],[295,1.0],[296,0.6],[297,1.2],[299,3.5333],[300,1.3333],[301,0.3333],[304,0.1333],[307,2.0667],[309,2.5333],[310,1.0],[311,2.1333],[313,0.4667],[316,2.3333],[319,0.9333],[320,7.8],[323,1.1333],[324,2.0],[325,1.5333],[326,2.0667],[329,0.1333],[330,2.8667],[332,0.9333],[333,0.4667],[341,0.1333],[342,1.0],[343,1.6],[344,0.6],[345,0.0667],[346,0.1333],[348,0.2667],[349,0.0667],[350,1.0667],[352,0.8667],[353,0.1333],[354,0.2],[355,2.8],[357,1.4667],[358,0.4],[359,0.9333],[362,2.8],[363,4.2667],[364,0.8667],[365,1.9333],[366,0.0667],[367,1.7333],[368,8.8667],[369,2.0],[370,2.4667],[373,0.8],[374,0.2667],[376,0.1333],[382,0.5333],[383,1.8667],[384,1.8],[385,1.6],[386,1.8667],[387,0.4667],[391,2.2],[395,2.6],[396,3.9333],[397,4.5333],[398,1.1333],[403,1.1333],[405,0.1333],[406,1.8],[407,2.0667],[411,1.2667],[412,0.9333],[415,0.2],[417,1.4],[420,2.6],[421,0.7333],[422,2.7333],[424,0.2],[425,1.1333],[426,0.8667],[427,0.8667],[428,6.0],[429,0.1333],[430,1.0667],[431,2.8],[432,1.0667],[433,1.0],[436,1.8],[437,0.0667],[438,0.8],[440,0.5333],[441,1.4667],[442,5.4],[443,2.2],[444,2.2],[445,3.4],[446,1.5333],[447,0.2667],[452,0.7333],[455,2.0667],[456,1.0],[457,0.8],[459,0.8],[460,0.1333],[461,2.0],[462,1.6],[467,1.4],[468,1.0],[469,1.0],[470,0.2667],[471,2.4],[473,0.0667],[474,3.2],[475,1.4],[476,0.4667],[478,0.3333],[479,3.6667],[480,3.2667],[481,0.3333],[482,1.4],[483,0.4],[486,0.0667],[487,0.2667],[488,3.5333],[490,6.8667],[491,2.6667],[492,0.2],[498,1.0],[499,0.7333],[500,0.2],[501,0.6],[502,0.9333],[504,0.4667],[506,0.6],[508,0.7333],[512,0.6667],[514,1.0],[517,1.2],[518,1.8],[520,0.6667],[521,0.8667],[522,1.5333],[525,0.1333],[527,0.6667],[528,0.2],[529,0.6667],[530,0.5333],[535,0.2667],[536,0.6],[539,0.6],[545,0.7333],[547,1.6],[548,0.4667],[550,1.0667],[553,0.3333],[554,0.6],[556,0.1333],[557,1.2667],[559,2.2],[561,1.6],[565,1.5333],[566,1.4667],[567,0.7333],[570,1.0],[571,0.1333],[572,0.1333],[574,0.1333],[575,0.4667],[578,1.4],[579,0.0667],[580,0.1333],[581,3.5333],[582,3.1333],[587,1.8],[588,2.0667],[589,1.4667],[590,0.2],[591,0.1333],[592,1.0],[594,1.0],[597,1.0],[598,1.6],[601,0.9333],[602,2.0],[606,0.4],[609,1.0],[610,0.9333],[614,0.7333],[617,1.4],[618,1.4667],[619,0.0667],[620,2.4667],[622,0.8],[623,0.0667],[628,1.4]],"ko":[[215,1.2667],[302,0.6667]],"hot":[105,138,143,151,160,161,165,166,199,202,212,225,244,279,299,339,356,357,363,373,385,391,394,397,411,426,432,433,436,445,467,468,476,478,487,488,512,538,553,561,569,581,582,612,620],"window":[101,629]},"clones":[["WT2","Wt",323,0],["WT2","Wt",20,0],["WT2","Wt",6,25],["WT2","Wt",5,4],["WT2","Wt",5,7],["WT2","Wt",4,0],["WT2","Wt",4,0],["WT2","Wt",4,5],["WT2","Wt",3,0],["WT2","Wt",3,20],["WT2","Wt",2,31],["WT2","Wt",2,9],["WT2","Wt",2,11],["WT2","Wt",2,11],["WT2","Wt",2,11],["WT2","Wt",2,14],["WT2","Wt",1,7],["WT2","Wt",1,0],["WT2","Wt",1,1],["WT2","Wt",1,18],["WT2","Wt",1,17],["WT2","Wt",1,16],["WT2","Wt",1,0],["WT2","Wt",1,0],["WT3","Wt",69,4],["WT3","Wt",30,6],["WT3","Wt",28,18],["WT3","Wt",19,25],["WT3","Wt",17,18],["WT3","Wt",15,8],["WT3","Wt",15,37],["WT3","Wt",13,14],["WT3","Wt",13,19],["WT3","Wt",12,10],["WT3","Wt",12,3],["WT3","Wt",11,14],["WT3","Wt",10,0],["WT3","Wt",9,0],["WT3","Wt",9,30],["WT3","Wt",8,21],["WT3","Wt",8,19],["WT3","Wt",8,11],["WT3","Wt",8,8],["WT3","Wt",6,20],["WT3","Wt",6,23],["WT3","Wt",5,0],["WT3","Wt",5,10],["WT3","Wt",5,5],["WT3","Wt",5,12],["WT3","Wt",4,12],["WT3","Wt",4,10],["WT3","Wt",1,14],["WT3","Wt",1,11],["WT3","Wt",1,12],["WT3","Wt",1,23],["WT3","Wt",1,0],["WT3","Wt",1,18],["WT4","Wt",30,0],["WT4","Wt",23,10],["WT4","Wt",23,5],["WT4","Wt",19,4],["WT4","Wt",19,0],["WT4","Wt",19,0],["WT4","Wt",19,0],["WT4","Wt",18,0],["WT4","Wt",18,2],["WT4","Wt",16,18],["WT4","Wt",14,0],["WT4","Wt",13,0],["WT4","Wt",13,2],["WT4","Wt",13,2],["WT4","Wt",11,0],["WT4","Wt",6,2],["KO1","Aicda-/-",64,0],["KO1","Aicda-/-",49,0],["KO1","Aicda-/-",19,0],["KO1","Aicda-/-",16,0],["KO1","Aicda-/-",15,0],["KO1","Aicda-/-",13,0],["KO1","Aicda-/-",12,0],["KO1","Aicda-/-",12,0],["KO1","Aicda-/-",10,0],["KO1","Aicda-/-",10,0],["KO1","Aicda-/-",6,0],["KO1","Aicda-/-",5,0],["KO1","Aicda-/-",5,0],["KO1","Aicda-/-",5,0],["KO1","Aicda-/-",1,0],["KO1","Aicda-/-",1,0],["KO2","Aicda-/-",16,0],["KO2","Aicda-/-",13,0],["KO2","Aicda-/-",12,0],["KO2","Aicda-/-",10,2],["KO2","Aicda-/-",10,0],["KO2","Aicda-/-",9,0],["KO2","Aicda-/-",6,0],["KO2","Aicda-/-",5,0],["KO2","Aicda-/-",5,0],["KO2","Aicda-/-",5,0],["KO2","Aicda-/-",5,0],["KO2","Aicda-/-",5,0],["KO2","Aicda-/-",5,0],["KO2","Aicda-/-",3,0],["KO2","Aicda-/-",2,0],["KO2","Aicda-/-",2,0],["KO3","Aicda-/-",31,0],["KO3","Aicda-/-",23,0],["KO3","Aicda-/-",15,0],["KO3","Aicda-/-",13,0],["KO3","Aicda-/-",10,0],["KO3","Aicda-/-",9,0],["KO3","Aicda-/-",9,0],["KO3","Aicda-/-",8,0],["KO3","Aicda-/-",7,0],["KO3","Aicda-/-",7,0],["KO3","Aicda-/-",7,0],["KO3","Aicda-/-",7,0],["KO3","Aicda-/-",6,0],["KO3","Aicda-/-",6,0],["KO3","Aicda-/-",6,0],["KO3","Aicda-/-",5,0],["KO3","Aicda-/-",5,0],["KO3","Aicda-/-",5,0],["KO3","Aicda-/-",5,0],["KO3","Aicda-/-",4,0],["KO3","Aicda-/-",4,0],["KO3","Aicda-/-",4,0],["KO3","Aicda-/-",3,0],["KO3","Aicda-/-",3,0],["KO3","Aicda-/-",3,0],["KO3","Aicda-/-",2,0],["KO3","Aicda-/-",2,0]],"largest":{"WT2":81.4,"WT3":19.2,"WT4":10.9,"KO1":26.3,"KO2":14.2,"KO3":14.8},"drawn":21}/*END*/;
  const W = 1400, H = 1376, o = [];
  const INK = '#1d2733', MID = '#4a5562', SOFT = '#6f7b88', OUTL = '#b4bec8', FAINT = '#cfd8df';
  const CASE = '#274b7e', CASEF = '#6d9ad0', CTRL = '#9a3439', CTRLF = '#d07a7e';
  const F = 'arial', TITLE = 22.5, HEAD = 20.5, SM = 17.5;

  const T = (text, x, y, e = {}) => { const t = Make.text(text, x, y, { fontSize: SM, color: INK, family: F, ...e }); o.push(t); return t; };
  const TC = (text, cx, y, e = {}) => { const t = T(text, cx, y, { align: 'center', ...e }); t.x = cx - t.w / 2; return t; };
  const R = (x, y, w, h, e = {}) => { const r = Make.rect(x, y, w, h, { stroke: 'none', strokeWidth: 0, radius: 0, ...e }); o.push(r); return r; };
  const E = (cx, cy, d, e = {}) => { const r = Make.ellipse(cx - d / 2, cy - d / 2, d, d, { stroke: 'none', strokeWidth: 0, ...e }); o.push(r); return r; };
  const card = (x, y, w, h, e = {}) => R(x, y, w, h, { fill: '#ffffff', stroke: OUTL, strokeWidth: 1.6, radius: 10, ...e });
  const LINK = (a, pa, b, pb, e = {}) => { const c = Make.connector({ id: a.id, port: pa }, { id: b.id, port: pb }, { color: MID, width: 2.2, headSize: 0.7, ...e }); o.push(c); return c; };
  const width = (s, e = {}) => Make.text(s, 0, 0, { fontSize: SM, family: F, ...e }).w;
  // a paragraph wrapped to maxW, with an optional bold lead-in on its first line; returns the y below it
  const para = (lead, text, x, y, maxW, e = {}, dy = 21.5) => {
    const leadW = lead ? width(lead + ' ', { bold: true }) : 0;
    const lines = [];
    let line = '';
    for (const word of text.split(' ')) {
      const next = line ? `${line} ${word}` : word;
      if (width(next, e) > maxW - (lines.length ? 0 : leadW) && line) { lines.push(line); line = word; } else line = next;
    }
    lines.push(line);
    if (lead) T(lead, x, y, { bold: true, color: e.leadColor || INK });
    lines.forEach((s, k) => T(s, x + (k ? 0 : leadW), y + k * dy, { color: MID, ...e }));
    return y + lines.length * dy;
  };
  const CODE = (data, body) => `D = json.loads(${JSON.stringify(JSON.stringify(data))})\n${body}`;

  // ───────── the run: inputs → SHMtrace → outputs ─────────
  const top = 12, th = 172;
  const inp = card(20, top, 320, th);
  T('You provide', 36, top + 12, { fontSize: HEAD, bold: true });
  // two small pileups: an AID-competent library (mutations) and an AID-null control (a stray error)
  const pile = (x, y, n, dots, col) => {
    for (let r = 0; r < n; r++) R(x, y + r * 11, 100, 5, { fill: FAINT, radius: 2.5 });
    dots.forEach(([f, r]) => E(x + f * 100, y + r * 11 + 2.5, 8, { fill: col }));
  };
  pile(36, top + 56, 4, [[0.2, 0], [0.6, 0], [0.4, 1], [0.8, 2], [0.3, 3]], CASE);
  pile(36, top + 112, 3, [[0.7, 1]], CTRL);
  T('reads per animal', 152, top + 48, { color: CASE, bold: true });
  T('any platform', 152, top + 69, { color: MID });
  T('AID-null controls', 152, top + 102, { color: CTRL, bold: true });
  T('two or more', 152, top + 123, { color: MID });
  T('+ a sample sheet', 36, top + 142, { color: MID });

  const core = card(380, top, 200, th, { fill: '#edf5f8', stroke: '#b3d0db' });
  o.push(Make.image(LOGO, 480 - 28, top + 16, 56, 56));
  TC('SHMtrace', 480, top + 80, { fontSize: HEAD, bold: true, color: '#1f6f8b' });
  TC('calls and tests', 480, top + 110, { color: MID });
  TC('against the controls', 480, top + 131, { color: MID });

  const outs = card(620, top, 760, th, { stroke: FAINT });
  T('You get', 636, top + 12, { fontSize: HEAD, bold: true });
  const tiles = [['report.html', 'every result on one page'], ['figures', 'PNG, SVG or PDF'],
    ['tables', 'Excel and CSV, every number'], ['project.json', 'reopens in the app']];
  const tw = (760 - 32 - 3 * 14) / 4;
  tiles.forEach(([name, text], k) => {
    const x = 636 + k * (tw + 14), y = top + 50, gx = x + 12, gy = y + 12;
    R(x, y, tw, 106, { fill: '#f4f6f8', radius: 8 });
    // a small glyph for each kind of output
    if (k === 0) { R(gx, gy, 34, 26, { fill: '#ffffff', stroke: MID, strokeWidth: 1.6, radius: 3 }); R(gx, gy, 34, 6, { fill: MID }); }
    if (k === 1) { [12, 20, 8, 16].forEach((h, j) => R(gx + 2 + j * 8, gy + 26 - h, 6, h, { fill: j % 2 ? CTRLF : CASEF })); R(gx, gy + 26, 34, 1.6, { fill: INK }); }
    if (k === 2) { R(gx, gy, 34, 26, { fill: '#ffffff', stroke: MID, strokeWidth: 1.6 }); [8.7, 17.3].forEach((d) => R(gx, gy + d, 34, 1.4, { fill: MID })); R(gx + 11, gy, 1.4, 26, { fill: MID }); }
    if (k === 3) T('{ }', gx + 2, gy, { fontSize: HEAD, bold: true, color: MID });
    T(name, gx + 44, gy + 4, { bold: true });
    para('', text, gx, y + 50, tw - 24, {}, 20.5);
  });
  LINK(inp, 'e', core, 'w');
  LINK(core, 'e', outs, 'w');

  // ───────── key for the example plots ─────────
  const ky = top + th + 20;
  const kt = 'Examples redrawn from one run (JH4 intron, nanopore, 500 reads per mouse):';
  T(kt, 20, ky, { color: MID });
  let kx = 20 + width(kt) + 16;
  E(kx + 5, ky + 10, 10, { fill: CASEF, stroke: CASE, strokeWidth: 1.6 }); kx += 16 + T('wild type', kx + 16, ky, { color: CASE }).w + 20;
  E(kx + 5, ky + 10, 10, { fill: CTRLF, stroke: CTRL, strokeWidth: 1.6 }); T('AID-null ({i|Aicda}−/−)', kx + 16, ky, { color: CTRL });

  // ───────── the six questions ─────────
  const gx0 = 20, cw = (W - 40 - 2 * 20) / 3, ch = 530, gy0 = ky + 42, rowGap = 20, ph = 240;
  const panel = (k, accent, title, file, plot, read, watch, also) => {
    const x = gx0 + (k % 3) * (cw + 20), y = gy0 + Math.floor(k / 3) * (ch + rowGap);
    card(x, y, cw, ch, { stroke: accent, strokeWidth: 1.6 });
    R(x, y, cw, 6, { fill: accent, radius: 3 });
    E(x + 32, y + 36, 34, { fill: accent });
    TC(String.fromCharCode(65 + k), x + 32, y + 25, { fontSize: HEAD, bold: true, color: '#ffffff' });
    T(title, x + 58, y + 16, { fontSize: HEAD, bold: true, color: accent });
    T(file, x + 58, y + 41, { color: SOFT, italic: true });
    o.push(PLOT({ name: `${String.fromCharCode(65 + k)}: ${file}`, x: x + 8, y: y + 64, w: cw - 16, h: ph, code: plot }));
    let ty = para('', read, x + 18, y + 64 + ph + 14, cw - 36, { color: INK });
    ty = para('Watch for:', watch, x + 18, ty + 10, cw - 36, { leadColor: accent });
    para('Also:', also.join(' · '), x + 18, ty + 12, cw - 36, { color: SOFT, leadColor: SOFT });
  };

  // 1. read accounting
  panel(0, MID, 'Did enough reads get through?', 'qc_waterfall', CODE({ qc: DATA.qc }, `
fig, ax = plt.subplots()
for name, ys in D["qc"].items():
    c = CASE if name.startswith("WT") else CTRL
    ax.plot(range(5), ys, color=c, linewidth=lw(0.9), marker="o", markersize=2.6, markeredgewidth=0, zorder=3)
ax.set_xticks(range(5), ["in", "quality,\\nlength", "on\\ntarget", "usable", "analysed"])
ax.set_xlim(-0.25, 4.25); ax.set_ylim(0, 105); ax.set_yticks([0, 50, 100])
ax.set_ylabel("% of reads left")
ax.grid(axis="y", color=FAINT, linewidth=lw(0.5)); ax.set_axisbelow(True)
`),
  'One line per library: the share of reads left after each filter. The last step is the depth you set, not a loss.',
  'a library far below the others.',
  ['coverage', 'quality_histogram', 'saturation']);

  // 2. the control background
  panel(1, CTRL, 'Are the calls real?', 'control_background', CODE({ bg: DATA.bg }, `
fig, ax = plt.subplots()
x, y = D["bg"]["pos"], D["bg"]["rate"]
ax.fill_between(x, 5e-5, y, color=CTRLF, alpha=0.75, linewidth=0, step="mid")
med = sorted(y)[len(y) // 2]
ax.axhline(med, color=INK, linestyle=(0, (3, 2)), linewidth=lw(0.7))
ax.text(x[-1] - 4, 1.8e-2, f"median {100 * med:.2g}%", ha="right", va="center", color=INK,
        bbox=dict(facecolor="white", edgecolor="none", pad=1.0))
ax.set_yscale("log"); ax.set_ylim(5e-5, 4e-2); ax.set_xlim(x[0], x[-1])
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{100 * v:g}%"))
ax.yaxis.set_minor_formatter(plt.NullFormatter())
ax.set_xlabel("Position in the amplicon (bp)"); ax.set_ylabel("Control error rate")
`),
  'The error rate the AID-null controls showed at each position. A mutation is called only if it clears this floor, with predicted false calls held under 5%.',
  'spikes: calls there need much more support.',
  ['volcano', 'single_read_excess', 'replicate_concordance']);

  // 3. how much: density per animal
  panel(2, '#1f6f8b', 'How much mutation?', 'density_by_group', CODE({ d: DATA.dens }, `
import math
fig, ax = plt.subplots()
for i, (g, e, f) in enumerate([("Wt", CASE, CASEF), ("Aicda-/-", CTRL, CTRLF)]):
    rows = [r for r in D["d"]["rows"] if r[1] == g]
    vals = [r[2] if r[2] > 0 else r[3] for r in rows]
    xs = [i - 0.12 + 0.12 * k for k in range(len(rows))]
    for xx, v, r in zip(xs, vals, rows):
        ax.scatter([xx], [v], s=area(5.2), facecolor=f if r[2] > 0 else "white", edgecolor=e, linewidth=lw(0.8), zorder=3)
    gm = math.exp(sum(math.log(v) for v in vals) / len(vals))
    ax.plot([i - 0.24, i + 0.24], [gm, gm], color=e, linewidth=lw(1.3), solid_capstyle="butt", zorder=4)
ax.set_yscale("log"); plain_log(ax.yaxis); ax.set_ylim(0.1, 2000)
ax.set_xticks([0, 1], ["Wt", "Aicda−/−"]); ax.get_xticklabels()[1].set_fontstyle("italic")
ax.set_xlim(-0.6, 1.6)
top = 400
ax.plot([0, 0, 1, 1], [top / 1.5, top, top, top / 1.5], color=INK, linewidth=lw(0.7))
ax.text(0.5, top * 1.3, f"P = {D['d']['p']:.4f},  {D['d']['fold']:.0f}×", ha="center", va="bottom", color=INK)
ax.set_ylabel("Unique mutations\\nper 10$^{5}$ bp")
`),
  'One dot per animal (open: none found); bars are geometric means. The test counts each animal once, however many reads it gave.',
  'frequency counts a clone\'s mutations in every read; unique density counts them once.',
  ['mutation_frequency_by_group', 'mutated_reads_by_group', 'mutation_load', 'combined_density', 'indels']);

  // 4. what kind: the spectrum
  const hot = Math.round(DATA.spec.gc_hot.reduce((a, b) => a + b) / DATA.spec.gc_hot.length);
  panel(3, '#6a4c93', 'What kind of mutation?', 'spectrum_by_group', CODE({ s: DATA.spec }, `
fig, ax = plt.subplots()
cl = D["s"]["classes"]
cols = ["#b9a8d3" if c[0] in "AT" else "#6a4c93" for c in cl]
ax.bar(range(12), D["s"]["group"], width=0.72, color=cols, zorder=2)
for a in D["s"]["animals"]:
    ax.scatter(range(12), a, s=area(2.2), facecolor="white", edgecolor=INK, linewidth=lw(0.5), zorder=3)
ax.set_xticks(range(12), [c.replace(">", "→") for c in cl], rotation=90)
ax.set_xlim(-0.6, 11.6); ax.set_ylim(0, 27); ax.set_yticks([0, 10, 20])
ax.set_ylabel("% of unique\\nmutations")
ax.text(2.5, 26.5, "from A or T", color="#8a73b3", va="top", ha="center")
ax.text(8.5, 26.5, "from C or G", color="#6a4c93", va="top", ha="center")
`),
  `The 12 base changes: bars for wild type, dots per mouse. AID leaves C→T and G→A at WRC/GYW motifs (${hot}% of G:C mutations here).`,
  'compare amplicons on the composition-corrected spectrum.',
  ['substitution_matrix', 'hotspot_targeting', 'spectrum_corrected', 'targeting_spectrum', 'mutability_hedgehog']);

  // 5. where: position along the amplicon
  panel(4, '#a3621b', 'Where along the amplicon?', 'positional_map', CODE({ p: DATA.pos }, `
fig, ax = plt.subplots()
wx, wy = zip(*D["p"]["wt"])
ax.bar(wx, wy, width=1.6, color=CASE, linewidth=0, zorder=2)
if D["p"]["ko"]:
    kx, ky = zip(*D["p"]["ko"])
    ax.bar(kx, [-v for v in ky], width=1.6, color=CTRL, linewidth=0, zorder=2)
ax.axhline(0, color=INK, linewidth=lw(0.6))
top = max(wy) * 1.25
ax.scatter(D["p"]["hot"], [top * 0.97] * len(D["p"]["hot"]), marker="|", s=area(4), color="#a3621b", linewidth=lw(0.6))
ax.text(D["p"]["window"][1], top * 0.88, "WRC/GYW hotspots", color="#a3621b", va="top", ha="right")
ax.text(D["p"]["window"][1], -top * 0.1, "AID-null", color=CTRL, ha="right", va="top")
ax.set_xlim(*D["p"]["window"]); ax.set_ylim(-top * 0.3, top)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(v):g}"))
ax.set_xlabel("Position in the amplicon (bp)"); ax.set_ylabel("% of reads\\nmutated")
`),
  'The share of reads mutated at each position: wild type up, AID-null down. Peaks on hotspot ticks are the signature of AID.',
  'a peak in the controls too: a germline variant or an artefact.',
  ['positional_by_sample', 'targeting_profile']);

  // 6. clones
  const big = Object.entries(DATA.largest).sort((a, b) => b[1] - a[1])[0];
  panel(5, '#3d5a80', 'Which clones carry them?', 'clone_map', CODE({ c: DATA.clones }, `
fig, ax = plt.subplots()
for name, g, n, m in D["c"]:
    wt = g == "Wt"
    ax.scatter([n], [m], s=area(3.6 if wt else 3.4), facecolor=CASEF if wt else "white", edgecolor=CASE if wt else CTRL,
               linewidth=lw(0.5), zorder=3 if wt else 2, alpha=0.85)
ax.set_xscale("log"); plain_log(ax.xaxis); ax.set_xlim(0.7, 600); ax.set_ylim(-2.5, 40)
big = max((r for r in D["c"] if r[1] == "Wt"), key=lambda r: r[2])
ax.annotate(f"{big[0]}: {big[2]} reads,\\nunmutated", (big[2], big[3] + 0.8), xytext=(big[2], 9),
            ha="center", va="bottom", color=CASE, arrowprops=dict(arrowstyle="-", color=CASE, linewidth=lw(0.6)))
ax.set_xlabel("Reads in the clone"); ax.set_ylabel("Mutations\\nin the clone")
`),
  `Clones from the V(D)J junction; each dot is one clone. In ${big[0]}, one unmutated clone held ${Math.round(big[1])}% of reads.`,
  'one clone repeating its mutations in many reads.',
  ['clonality', 'phasing_blocks']);

  // ───────── footer ─────────
  const fy = gy0 + 2 * ch + rowGap + 14;
  T('Other routes: a single sample from the command line writes spectrum, spatial and mutation_load; the legacy compare writes three more.', 20, fy, { color: SOFT });

  return { id: 'crm-figS21', name: 'Figure S21', width: W, height: H, background: '#ffffff', objects: o,
    alt: 'SHMtrace outputs: a project run takes reads and AID-null controls and writes a report, figures, tables and the full result; its 26 figures are grouped by six questions, each with an example from a J_H4 nanopore run and a note on how to read it.' };
}
