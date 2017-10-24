% GB Comment:
1.	100
2a. 100
2b. 100
2c. 100
3a 100
3b. 100 
3c. 90 An elaboration in the comparison of the genes would be appropriate.  	
Overall: 99

%HW3

%% Problem 1 - Smith-Waterman alignment
% Consider two sequences 'GTAATCC' and 'GTATCCG'

% Construct the scoring matrix for this with the parameters:
% match value = 2, mismatch value = -1, and gap penalty = -1. Use your
% solution to get the optimal alignment. If you prefer, it is acceptable to do this with
% pencil and paper, you can then take a snapshot of your solution and
% include it in your repository. 

matchval = 2; 
mismatchval = -1; 
ofdiag = ones(7) - eye(7); 
S = matchval*eye(7) + mismatchval*ofdiag;

seq1 = 'GTAATCC'; 
seq2 = 'GTATCCG'; 

[score, align, start] = swalign(seq1, seq2, 'Alphabet', 'nt', 'ScoringMatrix', S, 'Showscore', true, 'GapOpen', 1);   
align 
%% Problem 2 - using the NCBI databases and sequence alignments

% Erk proteins are critical signal transducers of MAP kinase signaling.
% Accessions numbers for ERK1 (also called MAPK3) and ERK2 (also called MAPK1) human mRNA are NM_002746 and
% NM_002745, respectively. 


% Part 1. Perform an alignment of the coding DNA sequences of ERK1 and
% ERK2. What fraction of base pairs in ERK1 can align to ERK2? 
Erk_accession = 'NM_002746'; 
Erk2_accession = 'NM_002745'; 
Erk_data = getgenbank(Erk_accession);
Erk2_data = getgenbank (Erk2_accession); 

Erk_seq = Erk_data.Sequence;
Erk2_seq = Erk2_data.Sequence; 

[score, align, start] = swalign(Erk_seq, Erk2_seq, 'Alphabet', 'nt',  'Showscore', true);  % I am using the default matrix and gap penality.  
align 

bp_aligned = length(strfind(align(2,:), '|')); % this gives the number of aligned bp. 
tot_bp = length(Erk_seq);

fract = bp_aligned/tot_bp
%this gives 0.5536

% Part2. Perform an alignment of the aminoacid sequences of ERK1 and ERK2.
% What fraction of amino acids align?

Erk_aa =Erk_data.CDS.translation
Erk2_aa = Erk2_data.CDS.translation
[score1, align1, start1] = swalign(Erk_aa, Erk2_aa, 'Alphabet', 'AA',  'Showscore', true);% I am using the default matrix and gap penality. 

aa_aligned = length(strfind(align1(2,:), '|')); % this gives the number of aligned aa. 
tot_aa = length(Erk_aa); %total number of Erk amino acids. 
frac_aa = aa_aligned/tot_aa 

% I get 0.8047 


% Part 3.  Use the NCBI tools to get mRNA sequences for the mouse genes ERK1 and
% ERK2 and align both the coding DNA sequences and protein sequences to the
% human versions. How similar are they? 

% searching online I got the accession number: 
% for ERK1 (MAPK3); NM_011952
% for ERK2: D10939

%%ERK HUMAN VS MOUSE
% DNA
Erk_m_data = getgenbank('NM_011952');
Erk2_m_data = getgenbank('D10939');

Erk_m_seq = Erk_m_data.Sequence;
Erk2_m_seq = Erk2_m_data.Sequence; 

[score2, align2, start2] = swalign(Erk_seq, Erk_m_seq, 'Alphabet', 'nt',  'Showscore', true);  % I am using the default matrix and gap penality.  
align2

bp_aligned2 = length(strfind(align2(2,:), '|')); % this gives the number of aligned bp. 
tot_bp2 = length(Erk_m_seq);
fract2 = bp_aligned2/tot_bp2
% I get 0.8521

%Amino acid
Erk_m_aa =Erk_m_data.CDS.translation
Erk2_m_aa = Erk2_m_data.CDS.translation
[score3, align3, start3] = swalign(Erk_aa, Erk_m_aa, 'Alphabet', 'AA',  'Showscore', true);% I am using the default matrix and gap penality. 

aa_aligned1 = length(strfind(align3(2,:), '|')); % this gives the number of aligned aa. 
tot_aa1 = length(Erk_m_aa); %total number of Erk amino acids. 
frac_aa1 = aa_aligned1/tot_aa1 
% I get 0.9658


%%ERK2 HUMAN VS MOUSE
% DNA
[score4, align4, start4] = swalign(Erk2_seq, Erk2_m_seq, 'Alphabet', 'nt',  'Showscore', true);  % I am using the default matrix and gap penality.  
align4

bp_aligned3 = length(strfind(align4(2,:), '|')); % this gives the number of aligned bp. 
tot_bp3 = length(Erk2_m_seq);
fract3 = bp_aligned3/tot_bp3
% I get ~0.6079

% Amino acid
[score5, align5, start5] = swalign(Erk2_aa, Erk2_m_aa, 'Alphabet', 'AA',  'Showscore', true);% I am using the default matrix and gap penality. 

aa_aligned2 = length(strfind(align5(2,:), '|')); % this gives the number of aligned aa. 
tot_aa2 = length(Erk2_m_aa); %total number of Erk amino acids. 
frac_aa2 = aa_aligned2/tot_aa2 

% I get 0.9916

%% Problem 3: using blast tools programatically

% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits. 


acession='NM_005359'; %acession needs to be a string
N = 4; % N needs to be a number 
accession_number = bigblast(acession, N); 
accession_number



% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 
[human non_human] = blastmatch (acession) %use accession definition above  

human 
non_human 
% if no human match, variable human stores the string 'no human entry'  

% Part 3. Choose at least one gene from the human genome and one gene from
% another organism and run your code from part 2 on the relevant accession
% numbers. Comment on the results. 

accession1 = 'AF228312' %human CD3zeta chain 
accession2 = 'LN515608' %GFP
[human1 non_human1] = blastmatch (accession1)
[human2 non_human2] = blastmatch (accession2)

% If i run my code for genes that are not from human, I get 'no human
% entry'. If the gene is from human origin, I always get some match from
% human. 


