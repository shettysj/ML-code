%Acquire data from text file and pad with zeros wherever data is missing
clear all


fid = fopen('MNISTnumImages5000.txt');
textLine = fgets(fid); % Read first line.
lineCounter = 1;
Data_orig=zeros(5000,784);
while ischar(textLine)
    
	% get into numbers array.
	numbers = sscanf(textLine, '%f ');
	% Put numbers into a cell array 

	for k = 1 : length(numbers)
        		Data_orig(lineCounter, k) = numbers(k);
	end
		
	% Read the next line.
    textLine = fgets(fid);
	lineCounter = lineCounter + 1;
end
fclose(fid);


% disp(ca);
