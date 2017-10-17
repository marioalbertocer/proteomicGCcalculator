#!/usr/bin/perl

use File::Basename;

print "\n\n\nPlease press 'c' plus enter for Chloroplasts, or 'm' for Mitochondria or 'b' for Bacteria: ";

my $answer = <STDIN>;
chomp $answer;

my $path = '/Users/marioalberto/Desktop/GCdataset';

if ($answer eq 'c') {
#	open (SALIDA, ">>$path/CGdinucleotide-chloro.txt");
	open (SALIDA, ">>$path/GC_withoutCGdinucleotides-chloro.txt");	
	opendir (DIR, "$path/mis cloroffn/");
} elsif ($answer eq 'm') {
	print "\n\n\nPlease choose one option \(1 or 2 \+ Enter\): \n\n1. Metozoan mitochondria\n2. Plant mitochondria\n\n\n";
	$answerM = <STDIN>;
	chomp $answerM;
	
	if ($answerM == 1) {
#		open (SALIDA, ">>$path/CGdinucleotide-Mitmet.txt");
		open (SALIDA, ">>$path/GC_withoutCGdinucleotides-Mitmet.txt");
		opendir (DIR, "$path/Mitocondria/mis mitmetffn/");
	} elsif ($answerM == 2) {
#		open (SALIDA, ">>$path/CGdinucleotide-Mitpla.txt");
		open (SALIDA, ">>$path/GC_withoutCGdinucleotides-Mitpla.txt");		
		opendir (DIR, "$path/Mitocondria/mis mitplaffn/");
	} else {
		die "run the script again and choose 1 or 2\n";
	}
} elsif ($answer eq 'b') {
	print "\n\n\nWrite the name of the taxa exactly as it appears in the genebank reports \(\.gbkfiles\) or \'all\' for all bacteria\:";
	$answer2 = <STDIN>;
	chomp $answer2;
#	open (SALIDA, ">>$path/CGdinucleotide_counter_report\-$answer2.txt");
	open (SALIDA, ">>$path/GC_withoutCGdinucleotides\-$answer2.txt");
	opendir (DIR, "$path/mis ffn/");
} else {
	die "run the script again and choose c, m or b\n";
}


my @dir = readdir DIR;
@dir = grep {/(\.ffn)|(\;)/} @dir;


if ($answer =~ /^(c|m)$/) {
	@dirffn = @dir;
} else {
	$_ = '' for my ($archivo); 
	@dirffn = '';
	foreach $archivo (@dir) {
		chomp $archivo;
		if ($answer2 eq 'all') {
			push (@dirffn, $archivo);
		} else {
			open(ARCHIVO, "$path/mis ffn/$archivo") or die "no se puede abrir $archivo:$!";
			@ffn = <ARCHIVO>;
			chomp @ffn;
			
			$gbk = $ffn[0] if $ffn[0] =~ /^>/;
			$gbk =~ s/^.*ref\|//g;
			$gbk =~ s/\.[0-9]{1,2}\|.*$//g;
			$gbk =~ s/\s//g;
			$gbk = "$gbk\.gbk";
			open (INFO, "$path/mis gbk/$gbk") or die "no se puede abrir \"$gbk\": $!";
			@info[0..15] = <INFO>;
			chomp @info;
			$taxonomia = join ('', @info[8..15]);
			$taxonomia =~ s/(^.*ORGANISM)|(REFERENCE.*$)//g;
			push (@dirffn, $archivo) if $taxonomia =~ $answer2; 
		}		
	}
}


if (@dirffn eq '') {die "Write a valid taxa\n\n"};
@dirffn = grep {/(\.ffn)|(\;)/} @dirffn;

$_ = '' for my ($archivoADN, $acceso);


#_________________________________________________________________________________________________________
	print "\n\nPor favor espere...\n\n";
#_________________________________________________________________________________________________________


foreach $archivoADN (@dirffn) {

	if ($answer eq 'c') {
		open(ARCHIVOADN, "$path/mis cloroffn/$archivoADN") or die "no se puede abrir $archivoADN:$!";
	} elsif ($answer eq 'm') {
		if ($answerM == 1) {
			open(ARCHIVOADN, "$path/mitocondria/mis mitmetffn/$archivoADN") or die "no se puede abrir $archivoADN:$!";
		} elsif ($answerM == 2) {
			open(ARCHIVOADN, "$path/mitocondria/mis mitplaffn/$archivoADN") or die "no se puede abrir $archivoADN:$!";
		}
	} else {
		open(ARCHIVOADN, "$path/mis ffn/$archivoADN") or die "no se puede abrir $archivoADN:$!";
	}

	print "$archivoADN\n\n";
	
	$count = 0;		
	$err = 0;
	my $gen;
#	@CGsyndinucl_total = ();
#	$CGsyndinucl_total = '';
#	@CGnonsyndinucl_total = ();
#	$CGnonsyndinucl_total = '';
	@GCsyn_total = ();
	@GCnonsyn_total = ();

	while (my $secuencia = <ARCHIVOADN>) {
    	 chomp $secuencia;                
		if ($secuencia =~ /^>/) {
			$acceso = $secuencia;
			$acceso =~ s/^.*ref\|//g;
			$acceso =~ s/\.[0-9]{1,2}\|.*$//g;
	 	       procesa_gen($gen)    
	 		       	if $gen;             
        		$gen = '';  
       			next; 
    		}
		$gen .= $secuencia;  
	}
	procesa_gen($gen) 
    	if $gen;

	sub procesa_gen {
    		my $gen = shift;
		$long_gen=length($gen);	
			
		if ((($answer eq 'b') & ($gen =~ /^[ACGT]TG/)) || (($answer eq 'm') || ($answer eq 'c'))) {
			$count++;	
			%SN_codigo = qw(AGA MGM AGG MGY CTA QTa CGA QGa CTG QTg CGG QGg GTA GTa GTG GTg TTA MTM TGA *** TTG MTY ATA ATM
					AAA AAM ATC ATQ AAC AAQ AGC AGQ AAG AAY ATT ATM AAT AAM AGT AGM ATG ATG ACA ACa ACC ACc ACG ACg
					ACT ACt CAA CAM CAC CAQ CAG CAY CAT CAM CCA CCa CTC CTc CCC CCc CGC CGc CCG CCg CTT CTt CCT CCt
					CGT CGt GAA GAM GAC GAQ GAG GAY GAT GAM GCA GCa GGA GGa GTC GTc GCC GCc GGC GGc GCG GCg GGG GGg
					GTT GTt GCT GCt GGT GGt TAA *** TTC TTQ TAC TAQ TGC TGQ TAG *** TTT TTM TAT TAM TGT TGM TGG TGG
					TCA TCa TCC TCc TCG TCg TCT TCt);
					
			if ($taxonomia =~ /(Entomoplasmatales)|(Mycoplasmataceae)|(Acholeplasmatales)|(Candidatus Hodgkinia cicadicola)|(Candidatus Zinderia insecticola CARI)|(Candidatus Nasuia deltocephalinicola str. NAS-ALF)/) {
				unless ($taxonomia =~ /laidlawii/) {
					$SN_codigo{TGA} = 'TGM'; 
					$SN_codigo{TGG} = 'TGY';
				}
			} elsif ($taxonomia =~ /(candidate division GN02)|(candidate division SR1)/) {
				$SN_codigo{TGA} = 'MGA'; 
				$SN_codigo{GGA} = 'YGa'; 
			} elsif ($answerM == 1) {
				$SN_codigo{AGA} = 'AGa'; 
				$SN_codigo{AGG} = 'AGg'; 
				$SN_codigo{AGC} = 'AGc'; 
				$SN_codigo{AGT} = 'AGt';
				$SN_codigo{ATG} = 'ATY'; 
			                         
				if ($archivoADN =~ /Metazoa.Cnidaria/) {
					$SN_codigo{TGA} = 'TGM'; 
					$SN_codigo{TGG} = 'TGY';
				} elsif ($archivoADN =~ /Metazoa;Platyhelminthes|Metazoa;Echinodermata/) {
					$SN_codigo{ATG} = 'ATG';
				} elsif ($archivoADN =~ /Chordata.Cephalochordata/) {
					$SN_codigo{AGA} = 'MGA'; 
					$SN_codigo{AGG} = 'AGY';
					$SN_codigo{AGC} = 'AGQ';
					$SN_codigo{GGA} = 'YGa';
					$SN_codigo{AGT} = 'AGM';
				} elsif ($archivoADN =~ /Chordata.Tunicata/) {
					$SN_codigo{AGA} = 'MGM'; 
					$SN_codigo{AGG} = 'MGY'; 
					$SN_codigo{AGC} = 'AGQ';
					$SN_codigo{GGA} = 'YGa';
					$SN_codigo{AGT} = 'AGM';
					$SN_codigo{GGG} = 'YGg';
				}
			}
				
			$SN = '';
			@SN = '';
			
			
			for ($i = 0; $i<$long_gen; $i += 3) {
				$codon = substr($gen, $i, 3);
				
				$SN .= $SN_codigo { $codon };
				
			}
			
			
			$GC = $SN;
#			$SN2 = $SN;
			
#			$SN =~ s/cG/i/g;
#			$SN =~ s/cg/i/g;
#			$SN =~ s/cY/i/g;
#			$SN =~ s/Cg/i/g;
#			$SN =~ s/Qg/i/g;
			
#			$SN2 =~ s/Cg/I/g;
#			$SN2 =~ s/CG/I/g;
#			$SN2 =~ s/CY/I/g;
#			$SN2 =~ s/cG/I/g;
#			$SN2 =~ s/QG/I/g;
			
			$GC =~ s/[cCQ][gGY]//g;		
			
#			@SN=split('', $SN);
#			@SN2=split('',$SN2); 
			@GC=split('', $GC);
			
			
			@GCsyn = grep {/[gc]/} @GC;
			@GCnonsyn = grep {/[GC]/} @GC;
#			@CGsyndinucl = grep {/i/} @SN;
#			@CGnonsyndinucl = grep {/I/} @SN2;
			
#			push @CGsyndinucl_total, @CGsyndinucl;
#			push @CGnonsyndinucl_total, @CGnonsyndinucl;
			push @GCsyn_total, @GCsyn;
			push @GCnonsyn_total, @GCnonsyn;
			
		} else {
			$err++;
		}
	}

#$CGsyndinucl_total = @CGsyndinucl_total;
#$CGnonsyndinucl_total = @CGnonsyndinucl_total;
$GCsyn_total = @GCsyn_total;
$GCnonsyn_total = @GCnonsyn_total;

#print "$CGsyndinucl_total\t";
#print "$CGnonsyndinucl_total\t";
print "$GCsyn_total\t";
print "$GCnonsyn_total\n";


# _______________________________Output________________________________________________________


#print SALIDA "\n$archivoADN\t$acceso\t$CGsyndinucl_total\t$CGnonsyndinucl_total\t$GCsyn_total\t$GCnonsyn_total\t";
print SALIDA "\n$archivoADN\t$acceso\t$GCsyn_total\t$GCnonsyn_total\t";

}