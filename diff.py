#!/usr/bin/env python

def comp_s(arr):
	s1 = arr[0]
	s2 = arr[1]
	out = ''
	k= 0
	
	if len(s1) != len(s2):
		raise Exception("DIFF SIZE")
	
	for i in range(len(s1)):
		c1 = s1[i]
		c2 = s2[i]
		
		if c1 == c2:
			out += ' '
		else:
			out += '|'
			k += 1
	
	print "- k",k
	print s1
	print out
	print s2
	print 

comp_s(['> chr21	42861433	+	chr21	42866505	-	n_discordant_reads=3/11,n_support=17/44	42	32	14	3	7	5	6	1.0	0.973354170381	chr21:42861433/42861434(+)->chr21:42866505/42866506(-):(discordant_mates:1,spanning_paired_1_t:1,spanning_paired_2:2,spanning_singleton_1:1)&chr21:42861433/42861434(+)->chr21:42870116/42870117(-):(spanning_paired_2:3,spanning_singleton_2_r:1)&chr21:42860320/42860321(+)->chr21:42861520/42861521(-):(discordant_mates:1,spanning_paired_2:1,spanning_singleton_1:1)&chr21:42860568/42860569(+)->chr21:42867210/42867211(-):(spanning_paired_1_t:2)&chr21:42859986/42859987(+)->chr21:42863878/42863879(-):(spanning_paired_1_t:1)&chr21:42860320/42860321(+)->chr21:42866505/42866506(-):(spanning_paired_2:1)&chr21:42866446/42866447(+)->chr21:42870110/42870111(-):(discordant_mates:1)',
        '< chr21	42861433	+	chr21	42866505	-	n_discordant_reads=3/11,n_support=17/44	42	32	14	3	7	5	6	1.0	0.973354170381	chr21:42861433/42861434(+)->chr21:42866505/42866506(-):(discordant_mates:1,spanning_paired_1_t:1,spanning_paired_2:2,spanning_singleton_1:1)&chr21:42861433/42861434(+)->chr21:42870116/42870117(-):(spanning_paired_2:3,spanning_singleton_2_r:1)&chr21:42860320/42860321(+)->chr21:42861520/42861521(-):(discordant_mates:1,spanning_paired_2:1,spanning_singleton_1:1)&chr21:42860568/42860569(+)->chr21:42867210/42867211(-):(spanning_paired_1_t:2)&chr21:42859986/42859987(+)->chr21:42863878/42863879(-):(spanning_paired_1_t:1)&chr21:42860320/42860321(+)->chr21:42866505/42866506(-):(spanning_paired_2:1)&chr21:42866446/42866447(+)->chr21:42870110/42870111(-):(discordant_mates:1)'])

comp_s(['> chr21	39875880	+	chr21	39876476	-	entropy=0.40563906223,n_discordant_reads=0/4,n_support=4/16	12	8	4	0	2	2	2	0.0	0.40563906223	chr21:39875880/39875881(+)-">chr21:39876476/39876477(-):(spanning_paired_2:3)&chr21:39874280/39874281(+)-">chr21:39876277/39876278(-):(spanning_paired_1_t:1)',
       '< chr21	39875880	+	chr21	39876476	-	entropy=0.40563906223,n_discordant_reads=0/4,n_support=4/16	12	8	4	0	2	2	2	0.0	0.40563906223	chr21:39875880/39875881(+)-">chr21:39876476/39876477(-):(spanning_paired_2:3)&chr21:39874280/39874281(+)-">chr21:39876277/39876278(-):(spanning_paired_1_t:1)'])

comp_s(['< chr21	39766042	+	chr21	39766803	-	entropy=0.0,n_discordant_reads=0/2,n_support=2/8	6	4	2	0	1	1	1	0.0	0.0	chr21:39766042/39766043(+)-">chr21:39766803/39766804(-):(spanning_paired_1_t:2)',
        '> chr21	39766042	+	chr21	39766803	-	entropy=0.0,n_discordant_reads=0/2,n_support=2/8	6	4	2	0	1	1	1	0.0	0.0	chr21:39766042/39766043(+)-">chr21:39766803/39766804(-):(spanning_paired_1_t:2)'])

comp_s(['< chr21	39797816	+	chr21	39803772	-	entropy=0.0,n_discordant_reads=0/2,n_support=2/8	6	4	2	0	1	1	1	0.0	0.0	chr21:39797816/39797817(+)-">chr21:39803772/39803773(-):(spanning_paired_1_t:2)',
       '> chr21	39797816	+	chr21	39803772	-	entropy=0.0,n_discordant_reads=0/2,n_support=2/8	6	4	2	0	1	1	1	0.0	0.0	chr21:39797816/39797817(+)-">chr21:39803772/39803773(-):(spanning_paired_1_t:2)'])

comp_s(['> chr21	39874520	+	chr21	39875907	-	n_discordant_reads=0/2,n_support=2/8	5	4	2	0	1	1	1	1.0	1.0	chr21:39874520/39874521(+)-">chr21:39875907/39875908(-):(spanning_paired_1_t:1,spanning_singleton_1:1)',
       '< chr21	39874520	+	chr21	39875907	-	n_discordant_reads=0/2,n_support=2/8	5	4	2	0	1	1	1	1.0	1.0	chr21:39874520/39874521(+)-">chr21:39875907/39875908(-):(spanning_paired_1_t:1,spanning_singleton_1:1)'])

comp_s(['> chr21	39844415	+	chr21	39845107	-	n_discordant_reads=1/4,n_support=2/16	3	1	1	1	2	2	2	0.0	1.0	chr21:39844415/39844416(+)-">chr21:39845107/39845108(-):(spanning_singleton_2_r:1)&chr21:39843961/39843962(+)-">chr21:39845291/39845292(-):(discordant_mates:1)',
       '< chr21	39844415	+	chr21	39845107	-	n_discordant_reads=1/4,n_support=2/16	3	1	1	1	2	2	2	0.0	1.0	chr21:39844415/39844416(+)-">chr21:39845107/39845108(-):(spanning_singleton_2_r:1)&chr21:39843961/39843962(+)-">chr21:39845291/39845292(-):(discordant_mates:1)'])

comp_s(['< chr21	39860702	-	chr21	42843890	+	entropy=0.0,n_support=2/8	2	0	0	2	1	1	1	0.0	0.0	chr21:39860702/39860703(-)-">chr21:42843890/42843891(+):(discordant_mates:2)',
       '> chr21	39860702	-	chr21	42843890	+	entropy=0.0,n_support=2/8	2	0	0	2	1	1	1	0.0	0.0	chr21:39860702/39860703(-)-">chr21:42843890/42843891(+):(discordant_mates:2)'])

comp_s(['< chr21	42839824	+	chr21	42840323	-	n_support=2/8	2	2	0	2	1	1	1	1.0	1.0	chr21:42839824/42839825(+)-">chr21:42840323/42840324(-):(discordant_mates:2)',
       '> chr21	42839824	+	chr21	42840323	-	n_support=2/8	2	2	0	2	1	1	1	1.0	1.0	chr21:42839824/42839825(+)-">chr21:42840323/42840324(-):(discordant_mates:2)'])

print "done"
