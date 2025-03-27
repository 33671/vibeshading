#ifdef GL_ES
precision mediump float;
#endif 
#define PI 3.1415926535
uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;
float sdCircle(vec2 position, float radius) {
    return length(position) - radius; 
}
float sawtooth(float x, float period) {
    return 2.0 * (mod(x, period) / period) - 1.0;
}

void main() {
    // 计算窗口的宽高比
    // rvberk
    float aspect = u_resolution.x / u_resolution.y;

    // 归一化坐标，并考虑宽高比
    vec2 position = gl_FragCoord.xy / u_resolution;
    // position.x *= aspect; // 调整 x 坐标以保持圆形
    position.x = (position.x - 0.5) * aspect + 0.5; // 调整 x 坐标以保持圆形居中

    vec2 mouse_position = u_mouse.xy / u_resolution;
    // mouse_position.x *= aspect; // 调整 x 坐标以保持圆形
    mouse_position.x = (mouse_position.x - 0.5) * aspect + 0.5; // 调整 x 坐标以保持圆形居中

    // 绘制圆形
    float circle_color = sdCircle(position - mouse_position,0.2);
    vec3 color = vec3(sin(circle_color * 200.0 - 2.0*PI*sawtooth(u_time,1.0)));
    gl_FragColor = vec4(color,1.0);
}