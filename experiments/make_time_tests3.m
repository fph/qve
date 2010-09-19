1;

resolution=80;

make_time_test3('latouche',linspace(0.7,1.4,resolution))
legend('NM','PI','PN');
matlab2tikz('E1-slides.tex','height','\\plotheight','width','\\plotwidth');

make_time_test3('latouche',linspace(0.82,0.9,resolution))
legend('NM','PI','PN');
matlab2tikz('E1zm-slides.tex','height','\\plotheight','width','\\plotwidth');

%make_time_test3('latouche2',linspace(0.2,0.95,resolution))
%legend('NM','PI','PN');
%matlab2tikz('E2-slides.tex','height','\\plotheight','width','\\plotwidth');

%make_time_test3('latouche2',linspace(0.32,0.38,resolution))
%legend('NM','PI','PN');
%matlab2tikz('E2zm-slides.tex','height','\\plotheight','width','\\plotwidth');

%make_time_test3('latouche2',linspace(0.83,0.85,resolution))
%legend('NM','PI','PN');
%matlab2tikz('E2zm2-slides.tex','height','\\plotheight','width','\\plotwidth');