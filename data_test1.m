function [ML0F,AP0F,tm,Fs] = data_test1()
ML0F=zeros(163,6000);
AP0F=zeros(163,6000);

[tm,signal,Fs,siginfo]=rdmat('data/BDS00001m');ML0F(1,:)=signal(:,7);AP0F(1,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00013m');ML0F(2,:)=signal(:,7);AP0F(2,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00028m');ML0F(3,:)=signal(:,7);AP0F(3,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00037m');ML0F(4,:)=signal(:,7);AP0F(4,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00055m');ML0F(5,:)=signal(:,7);AP0F(5,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00067m');ML0F(6,:)=signal(:,7);AP0F(6,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00073m');ML0F(7,:)=signal(:,7);AP0F(7,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00088m');ML0F(8,:)=signal(:,7);AP0F(8,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00097m');ML0F(9,:)=signal(:,7);AP0F(9,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00115m');ML0F(10,:)=signal(:,7);AP0F(10,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00127m');ML0F(11,:)=signal(:,7);AP0F(11,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00139m');ML0F(12,:)=signal(:,7);AP0F(12,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00145m');ML0F(13,:)=signal(:,7);AP0F(13,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00160m');ML0F(14,:)=signal(:,7);AP0F(14,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00172m');ML0F(15,:)=signal(:,7);AP0F(15,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00190m');ML0F(16,:)=signal(:,7);AP0F(16,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00202m');ML0F(17,:)=signal(:,7);AP0F(17,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00211m');ML0F(18,:)=signal(:,7);AP0F(18,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00223m');ML0F(19,:)=signal(:,7);AP0F(19,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00229m');ML0F(20,:)=signal(:,7);AP0F(20,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00250m');ML0F(21,:)=signal(:,7);AP0F(21,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00262m');ML0F(22,:)=signal(:,7);AP0F(22,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00271m');ML0F(23,:)=signal(:,7);AP0F(23,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00280m');ML0F(24,:)=signal(:,7);AP0F(24,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00292m');ML0F(25,:)=signal(:,7);AP0F(25,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00310m');ML0F(26,:)=signal(:,7);AP0F(26,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00319m');ML0F(27,:)=signal(:,7);AP0F(27,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00334m');ML0F(28,:)=signal(:,7);AP0F(28,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00346m');ML0F(29,:)=signal(:,7);AP0F(29,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00358m');ML0F(30,:)=signal(:,7);AP0F(30,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00370m');ML0F(31,:)=signal(:,7);AP0F(31,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00382m');ML0F(32,:)=signal(:,7);AP0F(32,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00391m');ML0F(33,:)=signal(:,7);AP0F(33,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00397m');ML0F(34,:)=signal(:,7);AP0F(34,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00415m');ML0F(35,:)=signal(:,7);AP0F(35,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00430m');ML0F(36,:)=signal(:,7);AP0F(36,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00439m');ML0F(37,:)=signal(:,7);AP0F(37,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00445m');ML0F(38,:)=signal(:,7);AP0F(38,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00457m');ML0F(39,:)=signal(:,7);AP0F(39,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00475m');ML0F(40,:)=signal(:,7);AP0F(40,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00481m');ML0F(41,:)=signal(:,7);AP0F(41,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00496m');ML0F(42,:)=signal(:,7);AP0F(42,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00511m');ML0F(43,:)=signal(:,7);AP0F(43,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00523m');ML0F(44,:)=signal(:,7);AP0F(44,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00538m');ML0F(45,:)=signal(:,7);AP0F(45,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00550m');ML0F(46,:)=signal(:,7);AP0F(46,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00553m');ML0F(47,:)=signal(:,7);AP0F(47,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00574m');ML0F(48,:)=signal(:,7);AP0F(48,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00583m');ML0F(49,:)=signal(:,7);AP0F(49,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00595m');ML0F(50,:)=signal(:,7);AP0F(50,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00607m');ML0F(51,:)=signal(:,7);AP0F(51,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00619m');ML0F(52,:)=signal(:,7);AP0F(52,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00628m');ML0F(53,:)=signal(:,7);AP0F(53,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00643m');ML0F(54,:)=signal(:,7);AP0F(54,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00649m');ML0F(55,:)=signal(:,7);AP0F(55,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00661m');ML0F(56,:)=signal(:,7);AP0F(56,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00679m');ML0F(57,:)=signal(:,7);AP0F(57,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00688m');ML0F(58,:)=signal(:,7);AP0F(58,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00697m');ML0F(59,:)=signal(:,7);AP0F(59,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00709m');ML0F(60,:)=signal(:,7);AP0F(60,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00730m');ML0F(61,:)=signal(:,7);AP0F(61,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00742m');ML0F(62,:)=signal(:,7);AP0F(62,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00748m');ML0F(63,:)=signal(:,7);AP0F(63,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00757m');ML0F(64,:)=signal(:,7);AP0F(64,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00769m');ML0F(65,:)=signal(:,7);AP0F(65,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00790m');ML0F(66,:)=signal(:,7);AP0F(66,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00796m');ML0F(67,:)=signal(:,7);AP0F(67,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00805m');ML0F(68,:)=signal(:,7);AP0F(68,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00826m');ML0F(69,:)=signal(:,7);AP0F(69,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00829m');ML0F(70,:)=signal(:,7);AP0F(70,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00850m');ML0F(71,:)=signal(:,7);AP0F(71,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00862m');ML0F(72,:)=signal(:,7);AP0F(72,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00865m');ML0F(73,:)=signal(:,7);AP0F(73,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00877m');ML0F(74,:)=signal(:,7);AP0F(74,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00892m');ML0F(75,:)=signal(:,7);AP0F(75,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00904m');ML0F(76,:)=signal(:,7);AP0F(76,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00919m');ML0F(77,:)=signal(:,7);AP0F(77,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00931m');ML0F(78,:)=signal(:,7);AP0F(78,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00943m');ML0F(79,:)=signal(:,7);AP0F(79,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00952m');ML0F(80,:)=signal(:,7);AP0F(80,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00969m');ML0F(81,:)=signal(:,7);AP0F(81,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00976m');ML0F(82,:)=signal(:,7);AP0F(82,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS00994m');ML0F(83,:)=signal(:,7);AP0F(83,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01006m');ML0F(84,:)=signal(:,7);AP0F(84,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01014m');ML0F(85,:)=signal(:,7);AP0F(85,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01021m');ML0F(86,:)=signal(:,7);AP0F(86,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01033m');ML0F(87,:)=signal(:,7);AP0F(87,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01048m');ML0F(88,:)=signal(:,7);AP0F(88,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01066m');ML0F(89,:)=signal(:,7);AP0F(89,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01072m');ML0F(90,:)=signal(:,7);AP0F(90,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01081m');ML0F(91,:)=signal(:,7);AP0F(91,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01093m');ML0F(92,:)=signal(:,7);AP0F(92,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01114m');ML0F(93,:)=signal(:,7);AP0F(93,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01117m');ML0F(94,:)=signal(:,7);AP0F(94,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01129m');ML0F(95,:)=signal(:,7);AP0F(95,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01144m');ML0F(96,:)=signal(:,7);AP0F(96,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01156m');ML0F(97,:)=signal(:,7);AP0F(97,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01168m');ML0F(98,:)=signal(:,7);AP0F(98,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01180m');ML0F(99,:)=signal(:,7);AP0F(99,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01192m');ML0F(100,:)=signal(:,7);AP0F(100,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01201m');ML0F(101,:)=signal(:,7);AP0F(101,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01219m');ML0F(102,:)=signal(:,7);AP0F(102,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01234m');ML0F(103,:)=signal(:,7);AP0F(103,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01246m');ML0F(104,:)=signal(:,7);AP0F(104,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01258m');ML0F(105,:)=signal(:,7);AP0F(105,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01267m');ML0F(106,:)=signal(:,7);AP0F(106,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01279m');ML0F(107,:)=signal(:,7);AP0F(107,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01294m');ML0F(108,:)=signal(:,7);AP0F(108,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01306m');ML0F(109,:)=signal(:,7);AP0F(109,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01315m');ML0F(110,:)=signal(:,7);AP0F(110,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01327m');ML0F(111,:)=signal(:,7);AP0F(111,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01339m');ML0F(112,:)=signal(:,7);AP0F(112,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01348m');ML0F(113,:)=signal(:,7);AP0F(113,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01357m');ML0F(114,:)=signal(:,7);AP0F(114,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01378m');ML0F(115,:)=signal(:,7);AP0F(115,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01390m');ML0F(116,:)=signal(:,7);AP0F(116,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01396m');ML0F(117,:)=signal(:,7);AP0F(117,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01414m');ML0F(118,:)=signal(:,7);AP0F(118,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01420m');ML0F(119,:)=signal(:,7);AP0F(119,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01438m');ML0F(120,:)=signal(:,7);AP0F(120,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01450m');ML0F(121,:)=signal(:,7);AP0F(121,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01453m');ML0F(122,:)=signal(:,7);AP0F(122,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01465m');ML0F(123,:)=signal(:,7);AP0F(123,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01483m');ML0F(124,:)=signal(:,7);AP0F(124,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01489m');ML0F(125,:)=signal(:,7);AP0F(125,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01504m');ML0F(126,:)=signal(:,7);AP0F(126,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01516m');ML0F(127,:)=signal(:,7);AP0F(127,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01534m');ML0F(128,:)=signal(:,7);AP0F(128,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01546m');ML0F(129,:)=signal(:,7);AP0F(129,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01558m');ML0F(130,:)=signal(:,7);AP0F(130,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01564m');ML0F(131,:)=signal(:,7);AP0F(131,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01582m');ML0F(132,:)=signal(:,7);AP0F(132,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01588m');ML0F(133,:)=signal(:,7);AP0F(133,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01600m');ML0F(134,:)=signal(:,7);AP0F(134,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01609m');ML0F(135,:)=signal(:,7);AP0F(135,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01626m');ML0F(136,:)=signal(:,7);AP0F(136,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01633m');ML0F(137,:)=signal(:,7);AP0F(137,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01645m');ML0F(138,:)=signal(:,7);AP0F(138,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01663m');ML0F(139,:)=signal(:,7);AP0F(139,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01669m');ML0F(140,:)=signal(:,7);AP0F(140,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01684m');ML0F(141,:)=signal(:,7);AP0F(141,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01693m');ML0F(142,:)=signal(:,7);AP0F(142,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01711m');ML0F(143,:)=signal(:,7);AP0F(143,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01726m');ML0F(144,:)=signal(:,7);AP0F(144,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01732m');ML0F(145,:)=signal(:,7);AP0F(145,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01741m');ML0F(146,:)=signal(:,7);AP0F(146,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01762m');ML0F(147,:)=signal(:,7);AP0F(147,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01771m');ML0F(148,:)=signal(:,7);AP0F(148,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01783m');ML0F(149,:)=signal(:,7);AP0F(149,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01789m');ML0F(150,:)=signal(:,7);AP0F(150,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01804m');ML0F(151,:)=signal(:,7);AP0F(151,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01822m');ML0F(152,:)=signal(:,7);AP0F(152,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01825m');ML0F(153,:)=signal(:,7);AP0F(153,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01837m');ML0F(154,:)=signal(:,7);AP0F(154,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01858m');ML0F(155,:)=signal(:,7);AP0F(155,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01867m');ML0F(156,:)=signal(:,7);AP0F(156,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01873m');ML0F(157,:)=signal(:,7);AP0F(157,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01888m');ML0F(158,:)=signal(:,7);AP0F(158,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01897m');ML0F(159,:)=signal(:,7);AP0F(159,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01915m');ML0F(160,:)=signal(:,7);AP0F(160,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01921m');ML0F(161,:)=signal(:,7);AP0F(161,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01939m');ML0F(162,:)=signal(:,7);AP0F(162,:)=signal(:,8);
[tm,signal,Fs,siginfo]=rdmat('data/BDS01951m');ML0F(163,:)=signal(:,7);AP0F(163,:)=signal(:,8);





end



