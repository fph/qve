1;
make_iter_test('latouche',0.5:0.01:2);
legend('NM','PI');
matlab2tikz('E1it.tex','height','\\plotheight','width','\\plotwidth');
make_iter_test('latouche2',0.3:0.01:0.9);
matlab2tikz('E2it.tex','height','\\plotheight','width','\\plotwidth');