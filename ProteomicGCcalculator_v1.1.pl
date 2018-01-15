#!/usr/bin/perl

#----------- PLEASE READ -------------

# In order to run this script you should download the file ProteomicGCcalculator.zip
# and unzip it. Open the file ProteomicGCcalculator_v1.1.pl and modify the path in the line
# my $path = '/Users/marioalberto/Desktop/ProteomicGCcalculator;'

# This script uses genomic databases of Bacteria, Archaea, plant chloroplasts and metazoan mitochondria. 
# These databases were downloaded from NCBI (ftp://ftp.ncbi.nlm.nih.gov/) in 2012. You can also 
# download more recent data, but then you need to re-build the paths like these:

# my $path = 'your_path/your_folder';

# For metazoan mitochondria files:
# your_path/your_folder/Mitochondria/my_mitmetffn
# your_path/your_folder/Mitochondria/my_mitmetgbk  

# For plant mitochondria files:
# your_path/your_folder/Mitochondria/my_mitplaffn 
# your_path/your_folder/Mitochondria/my_mitplagbk

# For chloroplasts files:
# your_path/your_folder/Chloroplasts/my_chloroffn
# your_path/your_folder/Chloroplasts/my_chlorogbk  

# For Archaea-Bacteria files:
# your_path/your_folder/Bacteria-Archaea/my_ffn
# your_path/your_folder/Bacteria-Archaea/my_gbk 

#---------------------------------------

use File::Basename;

# The program is going to ask some question in the terminal, you just need to answer 
# properly. Here, the script is filtering the data to process only the files you want.

print "\n\n\nPlease press 'c' plus enter for Chloroplasts, or 'm' for Mitochondria or",
		"'x' for Bacteria-Archaea: ";

my $answer = <STDIN>;
chomp $answer;

my $path = '/Users/marioceron/Desktop/ProteomicGCcalculator';

if ($answer eq 'c') {
	open (OUT, ">>$path/ChlorocGCreport.txt");
	opendir (DIR, "$path/Chloroplasts/my_chloroffn/");
} elsif ($answer eq 'm') {
	print "\n\n\nPlease choose one option \(1 or 2 \+ Enter\): \n\n1.", 
			"Metozoan mitochondria\n2. Plant mitochondria\n\n\n";
	$answerM = <STDIN>;
	chomp $answerM;
	
	if ($answerM == 1) {
		open (OUT, ">>$path/MitmetGCreport.txt");
		opendir (DIR, "$path/Mitochondria/my_mitmetffn/");
	} elsif ($answerM == 2) {
		open (OUT, ">>$path/MitplaGCreport.txt");
		opendir (DIR, "$path/Mitochondria/my_mitplaffn/");
	} else {
		die "run the script again and choose 1 or 2\n";
	}
} elsif ($answer eq 'x') {
	print "\n\n\nWrite the name of the taxa exactly as it appears in the genebank reports ",
	"\(\.gbkfiles\) or \'all\' for all bacteria\:";
	$answer2 = <STDIN>;
	chomp $answer2;
	open (OUT, ">>$path/ProteomicGCreport\-$answer2.txt");
	opendir (DIR, "$path/Bacteria-Archaea/my_ffn/");
} else {
	die "run the script again and choose c, m or b\n";
}

print "\n\n\nPlease choose one option \(1,2 or 3 \+ Enter\): \n\n1. GC calculations\n2. ",
		"Gene lenght calculations\n3. Both\n\n\n";
$answerC = <STDIN>;
chomp $answerC;

unless (($answerC == 1)||($answerC ==2)||($answerC == 3)) {
	die "run the script again and choose 1,2 or 3\n";
}

#_____________________________________________________________________________________
# In the next lines the script is going to use the .gbk files for getting the taxonomy.
# This part is very important for the step above
#_____________________________________________________________________________________

my @dir = readdir DIR;
@dir = grep {/(\.ffn)|(\;)/} @dir;


if ($answer =~ /^(c|m)$/) {
	@dirffn = @dir;
} else {
	$_ = '' for my ($file); 
	@dirffn = '';
	foreach $file (@dir) {
		chomp $file;
		if ($answer2 eq 'all') {
			push (@dirffn, $file);
		} else {
			open(FILE, "$path/Bacteria-Archaea/my_ffn/$file")
				or die "This file cannot be opened $file:$!";
			@ffn = <FILE>;
			chomp @ffn;
			
			$gbk = $ffn[0] if $ffn[0] =~ /^>/;
			$gbk =~ s/^.*ref\|//g;
			$gbk =~ s/\.[0-9]{1,2}\|.*$//g;
			$gbk =~ s/\s//g;
			$gbk = "$gbk\.gbk";
			open (INFO, "$path/Bacteria-Archaea/my_gbk/$gbk")
				or die "This file cannot be opened \"$gbk\": $!";
			@info[0..15] = <INFO>;
			chomp @info;
			$taxonomy = join ('', @info[8..15]);
			$taxonomy =~ s/(^.*ORGANISM)|(REFERENCE.*$)//g;
			push (@dirffn, $file) if $taxonomy =~ $answer2; 
		}		
	}
}

if (@dirffn eq '') {die "Write a valid taxa\n\n"};
@dirffn = grep {/(\.ffn)|(\;)/} @dirffn;

$_ = '' for my ($fileDNA, $acces);


#________________________________________________________________________________________
	print "\n\nPlease wait...\n\n";
#________________________________________________________________________________________


# ----------  Now the script is going open the files  ---------------


foreach $fileDNA (@dirffn) {
	
	$_ = '' for ($preTCB, $preTCBslrg, $preTC, $nonpreTC, $preTCslrg, $nonpreTCslrg);
	
	if ($answer eq 'c') {
		open(FILEDNA, "$path/Chloroplasts/my_chloroffn/$fileDNA")
			or die "This file cannot be opened $fileDNA:$!";
	} elsif ($answer eq 'm') {
		if ($answerM == 1) {
			open(FILEDNA, "$path/mitochondria/my_mitmetffn/$fileDNA")
				or die "This file cannot be opened $fileDNA:$!";
		} elsif ($answerM == 2) {
			open(FILEDNA, "$path/mitochondria/my_mitplaffn/$fileDNA")
				or die "This file cannot be opened $fileDNA:$!";
		}
	} else {
		open(FILEDNA, "$path/Bacteria-Archaea/my_ffn/$fileDNA")
			or die "This file cannot be opened $fileDNA:$!";
	}

	print "$fileDNA\t";
	
	$count = 0;		
	$err = 0;
	my $gen;
	@gen = '';
	@long_list = ();
	$A1 = 0; $T1 = 0; $C1 = 0; $G1 = 0;
	$A2 = 0; $T2 = 0; $C2 = 0; $G2 = 0;
	$A3 = 0; $T3 = 0; $C3 = 0; $G3 = 0;	
	$Count_A1 = 0;	$Count_C1 = 0;	$Count_G1 = 0;	$Count_T1 = 0;
	$Count_A2 = 0;	$Count_C2 = 0;	$Count_G2 = 0;	$Count_T2 = 0;
	$Count_A3 = 0;	$Count_C3 = 0;	$Count_G3 = 0;	$Count_T3 = 0;
	$N_A = 0;	$N_T = 0;	$N_C = 0;	$N_G = 0;	
	$N_a = 0;	$N_t = 0;	$N_c = 0;	$N_g = 0;	
	$Count_N_A = 0; $Count_N_C = 0; $Count_N_G = 0; $Count_N_T = 0;
	$Count_N_a = 0; $Count_N_c = 0; $Count_N_g = 0; $Count_N_t = 0;
	@F=(); @L=(); @I=(); @M=(); @V=(); @S=(); @P=(); @T=(); @A=(); @Y=(); @H=();
	@Q=(); @N=(); @K=(); @D=(); @E=(); @C=(); @W=(); @R=(); @G=(); @aa=();

	
# Here the script is reading each file. As these are .ffn files, their structure is like:
#----------- >tag \n gene_sequence >tag \n gene_sequence >tag \n gene_sequence ...-------
#----------- Then, the script contains a subrutin for defining a gene -------------------
	
	while (my $sequence = <FILEDNA>) {
    	 chomp $sequence;                
		if ($sequence =~ /^>/) {
			$acces = $sequence;
			$acces =~ s/^.*ref\|//g;
			$acces =~ s/\.[0-9]{1,2}\|.*$//g;
	 	       process_gen($gen)    
	 		       	if $gen;             
        		$gen = '';  
       			next; 
    		}
		$gen .= $sequence;  
	}
	process_gen($gen) 
    	if $gen;

	sub process_gen {
    		my $gen = shift;
		
		$long_gen=length($gen);	
		unless ($answerC == 1) {
			push (@long_list, $long_gen); 
		}


# -------- Translate each gene and be careful to reassign the different  -------------
#        codons in those particular cases that have different genetic codes
# --------- Then, apply the synonymous-non-synonymous code----------------------------
		
		unless ($answerC == 2) {  
			
			if ((($answer eq 'x') & ($gen =~ /^[ACGT]TG/)) || 
				(($answer eq 'm') || ($answer eq 'c'))) {
				
				$count++;	
				
				%SN_code = qw(AGA MGM AGG MGM CTA MTa CGA MGa CTG MTg CGG MGg 
						GTA GTa GTG GTg TTA MTM TGA *** TTG MTM ATA ATM
						AAA AAM ATC ATM AAC AAM AGC AGM AAG AAM ATT ATM 
						AAT AAM AGT AGM ATG ATG ACA ACa ACC ACc ACG ACg
						ACT ACt CAA CAM CAC CAM CAG CAM CAT CAM CCA CCa 
						CTC CTc CCC CCc CGC CGc CCG CCg CTT CTt CCT CCt
						CGT CGt GAA GAM GAC GAM GAG GAM GAT GAM GCA GCa 
						GGA GGa GTC GTc GCC GCc GGC GGc GCG GCg GGG GGg
						GTT GTt GCT GCt GGT GGt TAA *** TTC TTM TAC TAM 
						TGC TGM TAG *** TTT TTM TAT TAM TGT TGM TGG TGG
						TCA TCa TCC TCc TCG TCg TCT TCt);
						
				%Ge_code = qw(AGA R AGG R CTA L CGA R CTG L CGG R 
						GTA V GTG V TTA L TGA * TTG L ATA I
						AAA K ATC I AAC N AGC S AAG K ATT I 
						AAT N AGT S ATG M ACA T ACC T ACG T
						ACT T CAA Q CAC H CAG Q CAT H CCA P 
						CTC L CCC P CGC R CCG P CTT L CCT P
						CGT R GAA E GAC D GAG E GAT D GCA A 
						GGA G GTC V GCC A GGC G GCG A GGG G
						GTT V GCT A GGT G TAA * TTC F TAC Y 
						TGC C TAG * TTT F TAT Y TGT C TGG W
						TCA S TCC S TCG S TCT S);
						
				if ($taxonomy =~ /(Entomoplasmatales)|(Mycoplasmataceae)|(Acholeplasmatales)|(Candidatus Hodgkinia cicadicola)|(Candidatus Zinderia insecticola CARI)|(Candidatus Nasuia deltocephalinicola str. NAS-ALF)/) {
					unless ($taxonomy =~ /laidlawii/) {
					
						$SN_code{TGA} = 'TGM'; $Ge_code{TGA} = 'W';
						$SN_code{TGG} = 'TGM';
					}
				} elsif ($taxonomy =~ /(candidate division GN02)|(candidate division SR1)/) {
									
					$SN_code{TGA} = 'MGA'; $Ge_code{TGA} = 'G';
					$SN_code{GGA} = 'MGa';
					 
				} elsif ($answerM == 1) {
				
					$SN_code{AGA} = 'AGa';	$Ge_code{AGA} = 'S';
					$SN_code{AGG} = 'AGg';	$Ge_code{AGG} = 'S';
					$SN_code{AGC} = 'AGc'; 
					$SN_code{AGT} = 'AGt';
					$SN_code{ATG} = 'ATM'; 
											$Ge_code{ATA} = 'M';
					                         
					if ($fileDNA =~ /Metazoa.Cnidaria/) {
					
						$SN_code{TGA} = 'TGM'; $Ge_code{TGA} = 'W';
						$SN_code{TGG} = 'TGM';
						
					} elsif ($fileDNA =~ /Metazoa;Platyhelminthes|Metazoa;Echinodermata/) {
										
						$SN_code{ATG} = 'ATG';
						$Ge_code{ATA} = 'I';
						$Ge_code{AAA} = 'N';
									 
					} elsif ($fileDNA =~ /Chordata.Cephalochordata/) {
					
						$SN_code{AGA} = 'MGA'; $Ge_code{AGA} = 'G';
						$SN_code{AGG} = 'AGM';
						$SN_code{AGC} = 'AGM';
						$SN_code{GGA} = 'MGa';
						$SN_code{AGT} = 'AGM';
						
					} elsif ($fileDNA =~ /Chordata.Tunicata/) {
					
						$SN_code{AGA} = 'MGM'; $Ge_code{AGA} = 'G';
						$SN_code{AGG} = 'MGM'; $Ge_code{AGG} = 'G';
						$SN_code{AGC} = 'AGM';
						$SN_code{GGA} = 'MGa';
						$SN_code{AGT} = 'AGM';
						$SN_code{GGG} = 'MGg';
					}
				}
				
				$SN = '';
				@SN = '';

# -------- Here is the counting of the pre-termination and non-pretermination codons ----

				
				for ($i = 0; $i<$long_gen; $i += 3) {

					$codon = substr($gen, $i, 3);
					if ($codon =~ /^(TAT)|(TAC)|(TGT)|(TGC)|(TGG)|(CAA)|(CAG)|(AAA)|(AAG)|(GAA)|(GAG)$/) {
						$preTC ++;
					} elsif ($codon  =~ /^(TTC)|(TTT)|(ATT)|(ATC)|(ATG)|(ATA)|(CTT)|(CTG)|(CTA)|(CTC)|(CCT)|(CCC)|(CCA)|(CCG)|(ACT)|(ACC)|(ACA)|(ACG)|(GCT)|(GCC)|(GCA)|(GCG)|(CAT)|(CAC)|(AAT)|(AAC)|(GAT)|(GAC)$/) {
						$nonpreTC ++;
					} elsif ($codon  =~ /^(TCA)|(TCG)|(TTA)|(TTG)|(AGA)|(GGA)|(CGA)$/) {
						$preTCslrg ++;
					} elsif ($codon  =~ /^(TCT)|(TCC)|(AGT)|(AGC)|(CTA)|(CTG)|(CTC)|(CTT)|(CGT)|(CGC)|(CGG)|(AGG)|(GGT)|(GGC)|(GGG)$/) {
						$nonpreTCslrg ++;
					}
					
					$SN .= $SN_code { $codon };

#  --Saving the amino acids in lists ---> this step is important for the ENC calculation	
					
					if ($Ge_code{ $codon } eq 'F') {push @F, $codon}
					if ($Ge_code{ $codon } eq 'L') {push @L, $codon}
					if ($Ge_code{ $codon } eq 'I') {push @I, $codon}
					if ($Ge_code{ $codon } eq 'M') {push @M, $codon}
					if ($Ge_code{ $codon } eq 'V') {push @V, $codon}
					if ($Ge_code{ $codon } eq 'S') {push @S, $codon}
					if ($Ge_code{ $codon } eq 'P') {push @P, $codon}
					if ($Ge_code{ $codon } eq 'T') {push @T, $codon}
					if ($Ge_code{ $codon } eq 'A') {push @A, $codon}
					if ($Ge_code{ $codon } eq 'Y') {push @Y, $codon}
					if ($Ge_code{ $codon } eq 'H') {push @H, $codon}
					if ($Ge_code{ $codon } eq 'Q') {push @Q, $codon}
					if ($Ge_code{ $codon } eq 'N') {push @N, $codon}
					if ($Ge_code{ $codon } eq 'K') {push @K, $codon}
					if ($Ge_code{ $codon } eq 'D') {push @D, $codon}
					if ($Ge_code{ $codon } eq 'E') {push @E, $codon}
					if ($Ge_code{ $codon } eq 'C') {push @C, $codon}
					if ($Ge_code{ $codon } eq 'W') {push @W, $codon}
					if ($Ge_code{ $codon } eq 'R') {push @R, $codon}
					if ($Ge_code{ $codon } eq 'G') {push @G, $codon}
				}
				
				@SN=split('', $SN);

# ---- Here is the counting if the synonymous and non synonymous A, C, G, T ------------


				for($i = 0; $i < $long_gen; $i += 1) {

					$N_A = $SN[$i] =~ tr/A/D/;
				       	$Count_N_A = $Count_N_A + $N_A;	
					$N_T = $SN[$i] =~ tr/T/I/;
				       	$Count_N_T = $Count_N_T + $N_T;			
					$N_C = $SN[$i] =~ tr/C/O/;
				       	$Count_N_C = $Count_N_C + $N_C;		
					$N_G = $SN[$i] =~ tr/G/U/;
				       	$Count_N_G = $Count_N_G + $N_G;
					
					$N_a = $SN[$i] =~ tr/a/d/;
				       	$Count_N_a = $Count_N_a + $N_a;	
					$N_t = $SN[$i] =~ tr/t/i/;
				       	$Count_N_t = $Count_N_t + $N_t;		
					$N_c = $SN[$i] =~ tr/c/o/;
				       	$Count_N_c = $Count_N_c + $N_c;	
					$N_g = $SN[$i] =~ tr/g/u/;
				        $Count_N_g = $Count_N_g + $N_g;
				}
				
			} else {
				$err++;
			}
			

# ----- Here is the counting of the A, C, G, T in each position of the codons: 1,2,3 ----
			
			@gen=split('', $gen);
			
			for($i = 0; $i < $long_gen; $i=$i+3) {
				$A1 = $gen[$i] =~ tr/A/a/;
			       	$Count_A1 = $Count_A1 + $A1;	
				$C1 = $gen[$i] =~ tr/C/c/;
			       	$Count_C1 = $Count_C1 + $C1;			
				$G1 = $gen[$i] =~ tr/G/g/;
			       	$Count_G1 = $Count_G1 + $G1;		
				$T1 = $gen[$i] =~ tr/T/t/;
			       	$Count_T1 = $Count_T1 + $T1;
				
				$A2 = $gen[$i+1] =~ tr/A/a/;
			       	$Count_A2 = $Count_A2 + $A2;	
				$C2 = $gen[$i+1] =~ tr/C/c/;
			       	$Count_C2 = $Count_C2 + $C2;			
				$G2 = $gen[$i+1] =~ tr/G/g/;
			       	$Count_G2 = $Count_G2 + $G2;		
				$T2 = $gen[$i+1] =~ tr/T/t/;
			       	$Count_T2 = $Count_T2 + $T2;
				
				$A3 = $gen[$i+2] =~ tr/A/a/;
			       	$Count_A3 = $Count_A3 + $A3;	
				$C3 = $gen[$i+2] =~ tr/C/c/;
			       	$Count_C3 = $Count_C3 + $C3;			
				$G3 = $gen[$i+2] =~ tr/G/g/;
			       	$Count_G3 = $Count_G3 + $G3;		
				$T3 = $gen[$i+2] =~ tr/T/t/;
			       	$Count_T3 = $Count_T3 + $T3;
			}
		}
	}
	
	unless ($answerC == 2) {
		
# Here is a list of the lists of amino acids --> This is important for the ENC ---------		
		
		@aa = (\@F, \@L, \@I, \@M, \@V, \@S, \@P, \@T, \@A,
			   \@Y, \@H, \@Q, \@N, \@K, \@D, \@E, \@C, \@W, \@R, \@G);

# ----------- Here the ENC calculations start _______________________
		
		$_ = '' for ($Ns, $K2, $K3, $K4, $nj_2, $nj_3, $nj_4,
					 $nj_Fcfj_2, $nj_Fcfj_3, $nj_Fcfj_4);
		@X = ();
		
		foreach my $aa (@aa) {
		    
		    @X = @$aa;
		    
		    if ($X[0] =~ /[ATCG]{3}/) {
		        
		        div_families (%fam1, %fam2);
		        
		        my @fami = (\%fam1, \%fam2);
		        
		        foreach $fami (@fami) {  
    
		            my %fami = %$fami;  
		            my @m = keys %fami;
		            my $m = @m;
		            
		            if ($m > 0) { 
 
		                my $Fcfj = 0;
		                my @nj = values %fami;
		                my $nj = 0; $nj += $_ for @nj;
		             
		                foreach $ni (@nj) {
		                    $Fcfj += (($ni+1)/($nj + $m))**2;
		                }
		              
		                if ($m == 1)  {
		                    $Ns++;
		                } elsif ($m == 2) {
		                    $K2++;
		                    $nj_2 += $nj; 
		                    $nj_Fcfj_2 += $nj * $Fcfj;
		                } elsif ($m == 3) {
		                    $K3++;
		                    $nj_3 += $nj;
		                    $nj_Fcfj_3 += ($nj * $Fcfj);                    
		                } elsif ($m == 4) {
		                    $K4++;
		                    $nj_4 += $nj;
		                    $nj_Fcfj_4 += ($nj * $Fcfj);
		                }  
		            }
		        }
		    }
		}
		
		sub div_families {
		    
		    %ni_family = ();
		    %fam1 = ();
		    %fam2 = ();
		    my $control = '';
		    
		    foreach my $codon (@X) {
		        unless ($control =~ /$codon/) {
		            $control .= "$codon\t";
		            my @Fam_i = '';
		            @Fam_i = grep /$codon/, @X;
		            my $ni_Fi = @Fam_i;
		            $ni_family {$codon} = $ni_Fi;
		        }
		    }
			
		    my @family = keys %ni_family;
		    my $num_cod = @family;
			
		    if ($num_cod > 4) {

		        foreach my $cod_sinon (keys %ni_family) {
		            my $codonM = $cod_sinon;
		            $codonM =~ s/[ATCG]$//g;

		            if ($family[0] =~ /(^$codonM)/) {
		                $fam1 {$cod_sinon} = ($ni_family {$cod_sinon});

		            } else {

		                my @family_b = '';
		                if ($family_b[0] =~ /(^$codonM)|/) {
		                    $fam2 {$cod_sinon} = ($ni_family {$cod_sinon});

		                } else {
		                    die "This codon: $cod_sinon cannot be placed
								in any of the families";
		                }
		            }
		        }
		    } else {
		        foreach my $cod_sinon (keys %ni_family) {
		            $fam1 {$cod_sinon} = ($ni_family {$cod_sinon});
		        }
		    }
		}
		
		$Nc = $Ns;
		$Nc += (($K2*$nj_2)/$nj_Fcfj_2) if $nj_Fcfj_2 != 0;
		$Nc += (($K3*$nj_3)/$nj_Fcfj_3) if $nj_Fcfj_3 != 0;
		$Nc += (($K4*$nj_4)/$nj_Fcfj_4) if $nj_Fcfj_4 != 0;
		
# --------- The next lines are the results of the preTCB alculations --------------------
		
		$preTCB = (28/11)*($preTC/$nonpreTC);
        $preTCBslrg = ((15 * $preTCslrg)/(7 * $nonpreTCslrg));

# --- The results of genome size, proteome size Gcsyn, GCnonsyn, GC1, GC2, GC3 -------

		$proporcionbuenos = $count/($count + $err);	

		$proteomesize = $Count_A3 + $Count_C3 + $Count_G3 + $Count_T3;

		$GC_content_NS = (($Count_N_C + $Count_N_G) /
					($Count_N_C + $Count_N_G + $Count_N_A + $Count_N_T));

		$GC_content_S = (($Count_N_c + $Count_N_g) /
					($Count_N_c + $Count_N_g + $Count_N_a + $Count_N_t)); 
		
		$Counttotal = $Count_A1+$Count_A2+$Count_A3+$Count_C1+$Count_C2+$Count_C3

		+$Count_G1+$Count_G2+$Count_G3+$Count_T1+$Count_T2+$Count_T3;

		$CountGCtotal = $Count_C1+$Count_C2+$Count_C3+$Count_G1+$Count_G2+$Count_G3;

		$GC_content_total = $CountGCtotal/$Counttotal;
		
		$GC_content1 = (($Count_G1 + $Count_C1) /
					($Count_A1 + $Count_T1 + $Count_C1 + $Count_G1));

		$GC_content2 = (($Count_G2 + $Count_C2) /
					($Count_A2 + $Count_T2 + $Count_C2 + $Count_G2));

		$GC_content3 = (($Count_G3 + $Count_C3) /
					($Count_A3 + $Count_T3 + $Count_C3 + $Count_G3));
		
# -------------- Printing in the terminal ----------------------------
		
		print "Good genes: $count      ", "Bad genes: $err\n",
				"Proportion of good genes: $proporcionbuenos\n";
	
		print "\nProteome size: ", $proteomesize, "\n";  

		print "\nPosition 1: A: ", $Count_A1,"  C: ", $Count_C1,
				"G: ", $Count_G1,"  T: ", $Count_T1, "\n";

		print "GC-content1: ", $GC_content1, "\n";

		print "\nPosition 2: A: ", $Count_A2,"  C: ", $Count_C2,
				"G: ", $Count_G2,"  T: ", $Count_T2, "\n";

		print "GC-content2: ", $GC_content2, "\n";

		print "\nPosition 3: A: ", $Count_A3,"  C: ", $Count_C3,
				"G: ", $Count_G3,"  T: ", $Count_T3, "\n";

		print "GC-content3: ", $GC_content3, "\n\n";

		print "GC-content-total: ", $GC_content_total, "\n\n";

		print " - Nucleotides in synonymous positions - \n", "A: $Count_N_a     ", 
				"C: $Count_N_c     ", "G: $Count_N_g     ", "T: $Count_N_t\n\n";

		print " - Nucleotides in non-synonymous positions - \n", "A: $Count_N_A     ", 
				"C: $Count_N_C     ", "G: $Count_N_G     ", "T: $Count_N_T\n\n"; 
     
		print "GC-content in synonymous positions: ", $GC_content_S, "\n\n";

		print "GC-content in non-synonymous positions: ", $GC_content_NS, "\n\n";

		print "\-\-Pretermination codons\-\-\n";

		print "preTCB \= $preTCB \tpreTCBslrg = $preTCBslrg\n\n";
		
		print "\-\-Effective number of codons\-\-\n";

 		print "\nNs (=K1): $Ns\tK2: $K2\tK3: $K3\tK4: $K4\t";

		print "\n\nSum K2 nj x Fcfj\: $nj_Fcfj_2\n";

		print "Sum K3 nj x Fcfj\: $nj_Fcfj_3\n";

		print "Sum K4 nj x Fcfj\: $nj_Fcfj_4\n";

		print "Sum K2 nj\: $nj_2\n";

		print "Sum K3 nj\: $nj_3\n";

		print "Sum K4 nj\: $nj_4\n";

		print "ENC Nc\: $Nc\n\n";

	}
	
	unless ($answerC == 1) {
		
# --- These are the calculations of average gene size and the standard deviation --------
		
		$size_list = @long_list;
		
		my $sum = 0;
		$sum += $_ for @long_list;			
		$ave_long = $sum/$size_list;
		
		my $sum_dif = 0;
		$sum_dif += $_ for (my @difs = map {($_ - $ave_long)**2} @long_list);
		
		$d_stand = sqrt ($sum_dif/$size_list);
		
		print "\n\-\-Statistics\-\-\n";
		print "\nStandard Dev: $d_stand\tAver: $ave_long\n\n";
	}
	

	foreach $list_GCNS_gene (@list_GCNS_gene){
		print "$list_GCNS_gene\n";
	}

	print "\n\n\n";
# _______________________________Output_________________________________________________


	print OUT "\n$acces\t$fileDNA\t";
	
	if ($answerC != 2) {
		print OUT "$proteomesize\t$count\t$err\t$proporcionbuenos\t",
				"$GC_content_total\t$GC_content1\t$GC_content2\t",
				"$GC_content3\t$GC_content_S\t$GC_content_NS\t$preTCB\t",
				"$preTCBslrg\t$Nc\t";
	}
	
	if ($answerC != 1) {
		print OUT "$ave_long\t$d_stand";	
	}
}
