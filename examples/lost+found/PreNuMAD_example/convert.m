for i=1:length(blade.stations)
    tmp(i,:)=[blade.stations(i).spanlocation blade.stations(i).degreestwist blade.stations(i).chord blade.stations(i).percentthick];
end

save Convert tmp -ascii

blade.components(5).name
blade.components(5).cp