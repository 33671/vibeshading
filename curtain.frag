precision highp float;
uniform vec2 u_resolution;
uniform float u_time;
float sdSphere(vec3 position_3d,float radius)
{
    return length(position_3d)-radius;
}
float sdTorus(vec3 p,vec2 t)
{
    vec2 q=vec2(length(p.xz)-t.x,p.y);
    return length(q)-t.y;
}
float sdfBox( vec3 p, vec3 b ) {
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}
#define nsin(x) (sin(x) * 0.5 + 0.5)
// --- SDF Function for the Animated Curtain ---
float sdfCurtainInspired(vec3 p, float time) {
    // --- Parameters ---
    // Basic dimensions and position
    vec3 center = vec3(0.0, 0.0, 0.0); // Center of the curtain space
    vec2 sizeXY = vec2(100.0, 5.5);     // Width (X) and Height (Y)
    float thickness = .1;          // Thickness of the fabric

    // Wave parameters (inspired by 2D code, values may need tuning)
    float freq1 = 15.0; // Corresponds to 100.0 in 2D (adjust for 3D scale)
    float amp1 = 0.12;   // Corresponds to 0.075 in 2D (adjust for 3D scale)
    float speed1 = -1.0; // Time direction/speed for wave 1

    float freq2Base = 18.0; // Base frequency for wave 2 (adjust)
    float freq2DistScale = 15.0; // Scaling for distance effect on freq2 (adjust)
    float amp2 = 0.15;    // Corresponds to 0.1 in 2D (adjust)
    float speed2 = 1.2;  // Time direction/speed for wave 2

    // --- Calculations ---
    // 1. Transform point relative to curtain center
    vec3 rp = p - center; // Relative point

    // 2. Calculate vertical factor (0 at bottom, 1 at top) - for amplitude modulation
    // Map y from [-sizeXY.y/2, sizeXY.y/2] to [0, 1]
    float y_norm = clamp((rp.y + sizeXY.y * 0.5) / sizeXY.y, 0.0, 1.0);
    // Make folds deeper at the bottom (reverse of typical curtain, but matches aurora more)
    // Or make them shallower at the top (more like fabric):
    float verticalAmplitudeFactor = smoothstep(0.0, 0.8, 1.0 - y_norm); // Shallower folds higher up

    // 3. Calculate horizontal position normalized (e.g., -0.5 to 0.5) for distance calc
    float x_norm = rp.x / sizeXY.x; // Ranges approx -0.5 to 0.5 if centered

    // 4. Calculate wave displacement based on adapted 2D logic
    // Wave 1
    float wave1 = nsin(speed1 * time + rp.x * freq1) * (amp1 * verticalAmplitudeFactor);

    // Wave 2 (frequency depends on distance from center x)
    float distFromCenterX = abs(x_norm); // Distance from center normalized (0 to ~0.5)
    float currentFreq2 = freq2Base + distFromCenterX * freq2DistScale;
    float wave2 = nsin(speed2 * time + rp.x * currentFreq2) * (amp2 * verticalAmplitudeFactor);

    // Combine waves and potentially add offset (inspired by -0.5 in 2D)
    // The offset isn't strictly needed for shape, but can center the displacement
    float surfaceZ = wave1 + wave2 - (amp1 + amp2) * verticalAmplitudeFactor * 0.5; // Center displacement around 0

    // 5. Calculate distance to this undulating surface (primarily along Z)
    float distToSurface = rp.z - surfaceZ;

    // 6. Add thickness
    float dist = abs(distToSurface) - thickness * 0.5;

    // 7. Apply Bounds (intersect with a bounding box)
    // Z extent needs to accommodate max possible combined amplitude
    float maxWaveAmplitude = (amp1 + amp2); // Max possible peak amplitude
    vec3 boxHalfExtents = vec3(sizeXY.x * 0.5, sizeXY.y * 0.5, thickness * 0.5 + maxWaveAmplitude);
    float boxDist = sdfBox(rp, boxHalfExtents); // Use relative point rp

    // Intersect the undulating sheet distance with the bounding box distance
    dist = max(dist, boxDist);

    return dist;
}
float sawtooth(float x, float period) {
    return 2.0 * (mod(x, period) / period) - 1.0;
}
mat3 rotationMatrix(vec3 axis,float angle){
    axis=normalize(axis);
    float s=sin(angle);
    float c=cos(angle);
    float oc=1.-c;
    
    return mat3(oc*axis.x*axis.x+c,oc*axis.x*axis.y-axis.z*s,oc*axis.z*axis.x+axis.y*s,
        oc*axis.x*axis.y+axis.z*s,oc*axis.y*axis.y+c,oc*axis.y*axis.z-axis.x*s,
    oc*axis.z*axis.x-axis.y*s,oc*axis.y*axis.z+axis.x*s,oc*axis.z*axis.z+c);
}
vec3 rotatePoint(vec3 p, vec3 axis, float theta) {
    float cost = cos(theta);
    float sint = sin(theta);
    return p * cost + cross(axis, p) * sint + axis * dot(axis, p) * (1.0 - cost);
}
float smin(float a,float b,float k)
{
    k*=1.;
    float r=exp2(-a/k)+exp2(-b/k);
    return-k*log2(r);
}
vec3 sdfOrbit(vec3 p, float time) {
    float orbitRadius = 1.0;  
    float orbitSpeed = 1.0;   
    float spinSpeed = 2.0;    

   
    float thetaOrbit = time * orbitSpeed;
    vec3 orbitCenter = vec3(cos(thetaOrbit) * orbitRadius, 0.0, sin(thetaOrbit) * orbitRadius);

    vec3 localP = p - orbitCenter;

    float thetaSpin = time * spinSpeed;
    vec3 rotatedP = rotatePoint(localP, vec3(0.0, 1.0, 0.0), thetaSpin);
    return rotatedP;
}

vec3 colorTwinPeaksPattern(vec3 pos, float scale, float thickness, float chevronHeight) {
    // 使用x和z坐标作为平面坐标（忽略y坐标）
    vec2 p = pos.xz;
    
    // 应用缩放
    p *= scale;
    
    // 定义水平周期（一个完整V形的宽度）
    float periodX = 1.0;
    
    // 计算当前点在V形周期内的水平位置（-1到1）
    float vx = mod(p.x, periodX) / periodX * 2.0 - 1.0;
    
    // 计算V形带来的垂直偏移
    float vOffset = abs(vx) * chevronHeight;
    
    // 组合y坐标（实际上是z坐标）与V形偏移
    float combined_y = p.y + vOffset;
    
    // 定义垂直周期（两个波段的高度）
    float periodY = chevronHeight * 2.0;
    
    // 计算在双波段周期内的垂直位置（0到1）
    float tile_y = mod(combined_y, periodY) / periodY;
    
    // 计算带符号距离（与SDF函数相同）
    float signed_dist = mod(tile_y - thickness + 0.5, 1.0) - 0.5;
    
    // 根据距离决定颜色：黑色或黄色
    // 使用step函数创建硬边界，或者使用smoothstep创建抗锯齿边界
    float pattern = step(0.0, signed_dist);
    
    // 黑色和黄色
    vec3 black = vec3(0.0);
    vec3 yellow = vec3(0.85, 0.75, 0.65);;
    
    // 混合颜色
    vec3 color = mix(black, yellow, pattern);
    
    // 返回带不透明度的颜色
    return color;
}
float sdf_plane_xz(vec3 p)
{
    return p.y;
}
float sdf_scene(vec3 position_3d)
{
    mat3 rot = rotationMatrix(vec3(0, 1.0, 0.0), 3.14 / 2.0);
    float d0 = sdfCurtainInspired(rot*position_3d - vec3(-1.66,0,1.4),u_time*.4);
    float d = sdfCurtainInspired(position_3d - vec3(0,0,0),u_time*.4);
    float d1 = sdf_plane_xz(position_3d);

    return min(d0,min(d,d1));
    // return d;
}
vec3 calcNormal(vec3 p) {
    const float eps = 0.001;
    return normalize(vec3(
        sdf_scene(p + vec3(eps,0,0)) - sdf_scene(p - vec3(eps,0,0)),
        sdf_scene(p + vec3(0,eps,0)) - sdf_scene(p - vec3(0,eps,0)),
        sdf_scene(p + vec3(0,0,eps)) - sdf_scene(p - vec3(0,0,eps))
    ));
}
vec3 shading(vec3 hit_point, vec3 nom, vec3 rd) {
    
    // float rot_angle = sawtooth(u_time,3.0) * 2.0 * 3.1415; 
    // mat3 rot=rotationMatrix(vec3(0,.6,2),rot_angle);
    vec3 lightPos = vec3(-4, 6.0, 5.0);
    vec3 lightDir = normalize(lightPos - hit_point);
    
    // 漫反射
    float diff = max(dot(nom, lightDir), .0);
    vec3 zigzag_color = colorTwinPeaksPattern(hit_point, 2.5, .5, 0.35);
    vec3 diffuse = diff * vec3(0.4941, 0.0314, 0.0314);  // 橙色材质
    if (abs(hit_point.y) < 0.01) {
        diffuse = diff * zigzag_color;
    }
    // if (zigzag_color - vec3(0.0) > vec3(0.01)) {
    //     diff * vec3(0.4941, 0.0314, 0.0314);
    // }
    
    // 环境光
    vec3 ambient = vec3(0.1);
    
    // 镜面反射
    vec3 reflectDir = reflect(lightDir, nom);
    float spec = pow(max(dot(rd, reflectDir), 0.0), 50.0);
    vec3 specular = 0.5 * spec * vec3(1.0);
    
    return ambient + diffuse + specular;
}
vec3 background(vec3 rd) {
    float t = 0.5*(rd.y + 1.0);
    return mix(vec3(1.0), vec3(0.5, 0.7, 1.0), t);
}
vec3 camera_rd(vec2 position_xy,vec3 right_nom,vec3 up_nom,vec3 forward)
{
    float zoom=1.0;
    return normalize(position_xy.x*right_nom+position_xy.y*up_nom+zoom*forward);
}

vec3 set_camera(vec2 position_xy,vec3 cam_postion,vec3 object_position){
    vec3 forward=normalize(object_position-cam_postion);
    vec3 right=normalize(cross(vec3(0, 1, 0),forward));
    vec3 up=normalize(cross(forward,right));
    return camera_rd(position_xy,right,up,forward);
}
float ray_march(vec3 ro,vec3 rd_nom)
{
    vec3 next_point=ro;
    float total_d=0.;
    for(int i=0;i<1000;i++)
    {
        float new_d=sdf_scene(next_point);
        if(new_d<.001){
            return total_d;
        }
        next_point=next_point+(rd_nom*new_d*0.5);
        total_d+=new_d*0.5;
        if(total_d>100.){
            return 2e10;
        }
    }
    return 2e10;
}
vec3 hsv2rgb(vec3 c)
{
    vec4 K=vec4(1.,2./3.,1./3.,3.);
    vec3 p=abs(fract(c.xxx+K.xyz)*6.-K.www);
    return c.z*mix(K.xxx,clamp(p-K.xxx,0.,1.),c.y);
}

void main()
{
    float aspect=u_resolution.x/u_resolution.y;
    vec2 position=gl_FragCoord.xy/u_resolution;
    position.x=(position.x-.5)*aspect+.5;
    vec3 cam=vec3(-4,1,3) /1.0;
    vec3 rd=set_camera(position-vec2(.5),cam,vec3(0,0,0.5));
    float d = ray_march(cam,rd);
    vec3 hitpoint = cam + rd * d;
    vec3 color = shading(hitpoint,calcNormal(hitpoint),rd);
    if (d > 200.0){
        color = background(rd);
    }
    // vec3 rgb = pow(color, vec3(1.0/2.2)); // 近似伽马校正
    // vec3 color=hsv2rgb(vec3(d/.1));
    gl_FragColor=vec4(color,1);
}