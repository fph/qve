1;

resolution=40;

make_time_test('latouche',linspace(0.5,2,resolution))
legend('NM','PI');
matlab2tikz('E1tm.tex','height','\\plotheight','width','\\plotwidth');

make_time_test('latouche',linspace(0.82,1,resolution))
legend('NM','PI');
matlab2tikz('E1zm.tex','height','\\plotheight','width','\\plotwidth');

make_time_test('latouche2',linspace(0.2,0.95,resolution))
legend('NM','PI');
matlab2tikz('E2tm.tex','height','\\plotheight','width','\\plotwidth');

make_time_test('latouche2',linspace(0.32,0.38,resolution))
legend('NM','PI');
matlab2tikz('E2zm1.tex','height','\\plotheight','width','\\plotwidth');

make_time_test('latouche2',linspace(0.83,0.85,resolution))
legend('NM','PI');
matlab2tikz('E2zm2.tex','height','\\plotheight','width','\\plotwidth');