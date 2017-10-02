function [human non_human] = blastmatch (accession) 
 accession_number = bigblast(accession, 50);

human = 'No human entry' ;
for i = 1:length(accession_number) % accession_number is ranked from highest score to lowest. 
resulti = getgenbank(char(accession_number(i))); 
if string(resulti.Source) == "Homo sapiens (human)"
    human = resulti;
    disp(human);
    break  
end  
end 

if isstruct(human) == 0 
    disp('No human entry')
end 

for j = 1:length(accession_number)
    
    resultj = getgenbank(char(accession_number(j))); 
      
    if string(resultj.Source) ~= "Homo sapiens (human)"
        non_human =resultj; 
        disp(non_human); 
        break 
    else 
    end 
end 