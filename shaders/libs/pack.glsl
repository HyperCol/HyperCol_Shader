//----------------------------------------------------------------------------//

float PackUnorm2x4(vec2 xy) {
	return dot(floor(15.0 * xy + 0.5), vec2(1.0 / 255.0, 16.0 / 255.0));
}
float PackUnorm2x4(float x, float y) { return PackUnorm2x4(vec2(x, y)); }
vec2 UnpackUnorm2x4(float pack) {
	vec2 xy; xy.x = modf(pack * 255.0 / 16.0, xy.y);
	return xy * vec2(16.0 / 15.0, 1.0 / 15.0);
}
float UnpackUnorm2x4X(float pack) { return fract(pack * 255.0 / 16.0) * 16.0 / 15.0; }
float UnpackUnorm2x4Y(float pack) { return floor(pack * 255.0 / 16.0) / 15.0; }

float PackUnorm2x8(vec2 xy) {
	return dot(floor(255.0 * xy + 0.5), vec2(1.0 / 65535.0, 256.0 / 65535.0));
}
float PackUnorm2x8(float x, float y) { return PackUnorm2x8(vec2(x, y)); }
vec2 UnpackUnorm2x8(float pack) {
	vec2 xy; xy.x = modf(pack * 65535.0 / 256.0, xy.y);
	return xy * vec2(256.0 / 255.0, 1.0 / 255.0);
}
float UnpackUnorm2x8X(float pack) { return fract(pack * 65535.0 / 256.0) * 256.0 / 255.0; }
float UnpackUnorm2x8Y(float pack) { return floor(pack * 65535.0 / 256.0) / 255.0; }

float PackSnorm2x8(vec2 xy) {
	vec2 xy2 = floor(mix(128.5 + 127.0 * xy, 127.5 + 127.0 * xy, greaterThan(xy, vec2(0.0))));
	return dot(xy2, vec2(1.0 / 65535.0, 256.0 / 65535.0));
}
float PackSnorm2x8(float x, float y) { return PackSnorm2x8(vec2(x, y)); }
vec2 UnpackSnorm2x8(float pack) {
	vec2 xy; xy.x = modf(pack * 65535.0 / 256.0, xy.y);
	return xy * vec2(256.0 / 127.0, 1.0 / 127.0) - mix(vec2(1.0), vec2(128.0 / 127.0), greaterThan(xy, vec2(127.5 / 256.0, 127.5)));
}
float UnpackSnorm2x8X(float pack) { return UnpackSnorm2x8(pack).x; }
float UnpackSnorm2x8Y(float pack) { return UnpackSnorm2x8(pack).y; }

//----------------------------------------------------------------------------//