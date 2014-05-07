function lg = bin2m(filename)

	fid = fopen(filename, 'r');
	dim = fread(fid, 1, 'int');
	
	lg = fread(fid, dim, 'int').';

	fclose(fid);
