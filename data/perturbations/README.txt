# parse_perturbation_file

Rules:
Accepts 1 file at a time using the option --pert_file_path.
File can have 1 or more perturbations.
Each perturbation occupies 1 row.
Tab-delimit different categories (e.g., lambda and 0 should be lambda\tab0).
Separate by commas if items belong to the same category (e.g., lambda and mu should be lambda,mu).

Headings can be in any order, as long as the values below match the headings.

Headings can be in upper/lower/mixed cases.

Headings' spelling can be “params”, “parameter”, “p” etc. (within reasonable variations)
Same goes for “values”, “update_mode”, and “axes”.

If a heading’s name cannot be recognized, an exception will be raised.
For the perturbations files, rows under “parameters” and “values” now accept
one or more inputs. If the number of parameters does not match the number of
values, then an exception will be raised.

If any required headings are missing, an exception will be raised.


Algorithm:
Each line of entry is put into a dictionary.
Each dictionary is put into a perturbations list.
The perturbations list is put into a list of treatments.


EXAMPLE Format:
parameters	values	update_mode	axes
lambda	0	replace	x,y,z
mu	0.8	add	x


EXAMPLE Format:
parameters	values	update_mode	axes
lambda,mu	0.016,-0.08	replace	x
