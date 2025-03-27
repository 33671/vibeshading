precision highp float;
uniform vec2 u_resolution;
uniform float u_time;

// 定义立方体距离函数
float sdBox(vec3 p, vec3 b) {
    vec3 q = abs(p) - b;
    return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

// 光线行进函数
float rayMarch(vec3 ro, vec3 rd) {
    float d = 0.0;
    for (int i = 0; i < 100; i++) {
        vec3 p = ro + rd * d;
        float dist = sdBox(p, vec3(0.5));
        if (dist < 0.001 || d > 100.0) break;
        d += dist;
    }
    return d;
}

// 计算法线
vec3 calcNormal(vec3 p) {
    const float h = 0.0001;
    vec2 k = vec2(1.0, -1.0);
    return normalize(k.xyy * sdBox(p + k.xyy * h, vec3(0.5)) + 
                     k.yyx * sdBox(p + k.yyx * h, vec3(0.5)) + 
                     k.yxy * sdBox(p + k.yxy * h, vec3(0.5)) + 
                     k.xxx * sdBox(p + k.xxx * h, vec3(0.5)));
}

// 主函数
vec4 mainImage(in vec2 fragCoord) {
    // 标准化坐标
    vec2 uv = (fragCoord - 0.5 * u_resolution.xy) / u_resolution.y;
    
    // 相机设置
    vec3 ro = vec3(0.0, 0.0, 2.0);
    vec3 rd = normalize(vec3(uv, -1.0));
    
    // 旋转立方体
    float angle = u_time;
    mat2 rot = mat2(2.*cos(angle), -sin(angle), sin(angle),2.0* cos(angle));
    rd.xy = rot * rd.xy;
    ro.xy = rot * ro.xy;
    
    // 光线行进
    float d = rayMarch(ro, rd);
    
    // 着色
    vec3 col = vec3(0.0);
    if (d < 100.0) {
        vec3 p = ro + rd * d;
        vec3 normal = calcNormal(p);
        
        // 简单光照
        vec3 lightPos = vec3(1.0, .5, 1.0);
        vec3 lightDir = normalize(lightPos - p);
        float diff = max(dot(normal, lightDir), 0.0);
        col = vec3(0.5, 0.7, 1.0) * diff;
        
        // 添加一些环境光
        col += vec3(0.1, 0.5, 0.2);
    }
    
   return vec4(col, 1.0);
}

void main() {
    gl_FragColor = mainImage(gl_FragCoord.xy);
}