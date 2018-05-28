%This sript is made by ROB18_Gr465 AAU 

%This scipt run the Playrec function which gather the data from the 3
%microphones[1 2 3] and with a sampling freqency of 96 kHz.
%When this  script is run the program asks the user for what driver to use
%-> pick the audiobox 14 channels on for both play_device and rec_device
%After the playrec script is done the multilateration script will run and
%finds all TDOAs from all the recordings and calculate the hyperbolas.
%futhermore is all intersection found between all hyperbolas to find the
%object at hand. the place where 3 intersection in the same place is
%considered a object.
line = Play_rec_new3(select_play_device(),select_rec_device(),1,[1 2 3],96000);
New_test(line)
