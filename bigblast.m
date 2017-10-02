function accession_number = bigblast (accession, N) 
gb_data = getgenbank(accession);
seq = gb_data.Sequence(1:200);
[requestID, requestTime] = blastncbi(seq, 'blastn');
blast_data = getblast(requestID, 'WaitTime', requestTime); 

accession_number = [ ]; 
for i = 1:N 
    name = blast_data.Hits(i).Name; 
    indices = strfind(name, '|'); 
    last = indices(end) -1;
    first = indices(end-1)+1; 
    accession_number{i,1} = name(first:last) ;
   
    
end 