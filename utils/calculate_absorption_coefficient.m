function [alpha, area] = calculate_absorption_coefficient(nRooms,room_dims, T60)
%% Calculate absorption coefficient of a room according to Sabine's formula
% INPUTS:
% nRooms - number of rooms to calculate for
% room_dims - nRoomsx3 room dimensions in metre
% T60 - nRoomsx1 T60 of each room
% OUTPUTS:
% alpha - absorption coefficient
% area - surface area of each room

alpha = zeros(nRooms,1);
area =  zeros(nRooms,1);
for k = 1:nRooms 
    volume = prod(room_dims(k,:));
    area(k) = 2*(room_dims(k,1)* room_dims(k,2) + room_dims(k,2)*room_dims(k,3) + room_dims(k,3)*room_dims(k,1));
    alpha(k) = 0.161 * volume/(T60(k)*area(k));
end