# SimpleFOC  v2.0.1 库解读 

## 一、库文件分析

## 1.common 基础文件

文件里是基本的foc底层算法，pid，低通滤波，三角函数等

### 1）defaults.h 

配置默认的参数值，可通过程序初始化设定进行改动

```C++
//电源电压
DEF_POWER_SUPPLY 12.0 
    
//速度环PID控制器参数（速度环通常只用到PI）
DEF_PID_VEL_P 0.5 
DEF_PID_VEL_I 10.0 
DEF_PID_VEL_D 0.0 
DEF_PID_VEL_U_RAMP 1000.0 //默认PID控制器电压斜坡值
    
//位置环P控制器参数（位置环这里只用到P）
DEF_P_ANGLE_P 20.0 
DEF_VEL_LIM 20.0 //默认角速度限制
    
//索引搜索目标速度 
DEF_INDEX_SEARCH_TARGET_VELOCITY 1.0 //默认索引搜索速度

//校准电压
DEF_VOLTAGE_SENSOR_ALIGN 6.0 //传感器和电机校准调零的默认电压 

//低通滤波速度
DEF_VEL_FILTER_Tf 0.005 //默认低通滤波速度时间常数
```

### 2）foc_utils

头文件定义：

```c
// 函数符号
_sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )
_round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
_constrain(amt,low,high) ((amt)<(low)?(low):((amt)>(high)?(high):(amt)))//约束

// 基本参数定义
#define _2_SQRT3 1.15470053838
#define _SQRT3 1.73205080757
#define _1_SQRT3 0.57735026919
#define _SQRT3_2 0.86602540378
#define _SQRT2 1.41421356237
#define _120_D2R 2.09439510239
#define _PI 3.14159265359
#define _PI_2 1.57079632679
#define _PI_3 1.0471975512
#define _2PI 6.28318530718
#define _3PI_2 4.71238898038

#define NOT_SET -12345.0
```

**<u>角度相关计算函数：</u>**

> 使用int数组代替float数组进行运算，2x存储保存（int 2Byte  float 4Byte）
>
> sin*10000

```c++
int sine_array[200] = {0,79,158,237,316,395,473,552,631,710,789,867,946,1024,1103,1181,1260,1338,1416,1494,1572,1650,1728,1806,1883,1961,2038,2115,2192,2269,2346,2423,2499,2575,2652,2728,2804,2879,2955,3030,3105,3180,3255,3329,3404,3478,3552,3625,3699,3772,3845,3918,3990,4063,4135,4206,4278,4349,4420,4491,4561,4631,4701,4770,4840,4909,4977,5046,5113,5181,5249,5316,5382,5449,5515,5580,5646,5711,5775,5839,5903,5967,6030,6093,6155,6217,6279,6340,6401,6461,6521,6581,6640,6699,6758,6815,6873,6930,6987,7043,7099,7154,7209,7264,7318,7371,7424,7477,7529,7581,7632,7683,7733,7783,7832,7881,7930,7977,8025,8072,8118,8164,8209,8254,8298,8342,8385,8428,8470,8512,8553,8594,8634,8673,8712,8751,8789,8826,8863,8899,8935,8970,9005,9039,9072,9105,9138,9169,9201,9231,9261,9291,9320,9348,9376,9403,9429,9455,9481,9506,9530,9554,9577,9599,9621,9642,9663,9683,9702,9721,9739,9757,9774,9790,9806,9821,9836,9850,9863,9876,9888,9899,9910,9920,9930,9939,9947,9955,9962,9969,9975,9980,9985,9989,9992,9995,9997,9999,10000,10000};
```

> 正弦函数（用定长数组逼近正弦函数），输入0-2π角度

```C++
float _sin(float a){
  if(a < _PI_2){
    return 0.0001*sine_array[_round(126.6873* a)];       //整数数组优化
  }else if(a < _PI){
    return 0.0001*sine_array[398 - _round(126.6873*a)];  
  }else if(a < _3PI_2){
    return -0.0001*sine_array[-398 + _round(126.6873*a)]; 
  } else {
    return -0.0001*sine_array[796 - _round(126.6873*a)];  
  }
}
```

> 余弦函数（用定长数组逼近正弦函数），输入0-2π角度

```c
float _cos(float a){
  float a_sin = a + _PI_2;
  a_sin = a_sin > _2PI ? a_sin - _2PI : a_sin;
  return _sin(a_sin);
}
```

> 弧度规格化成0-2π

```c
float _normalizeAngle(float angle){
  float a = fmod(angle, _2PI);
  return a >= 0 ? a : (a + _2PI);
}
```

> 电角度计算（shaft_angle 轴角， pole_pairs 极对数）

```c
float _electricalAngle(float shaft_angle, int pole_pairs) {
  return (shaft_angle * pole_pairs);
}
```



### 3）lowpass_filter

<u>**低通滤波**</u>

> 头文件类定义：

```c++
class LowPassFilter
{
public:
    /**
     * @param Tf - Low pass filter time constant
     */
    LowPassFilter(float Tf);
    ~LowPassFilter() = default;

    float operator() (float x);
    float Tf; //!< Low pass filter time constant

protected:
    unsigned long timestamp_prev;  //!< Last execution timestamp
    float y_prev; //!< filtered value in previous execution step 
};
```

> 函数：

```c++
LowPassFilter::LowPassFilter(float time_constant)
    : Tf(time_constant)
    , y_prev(0.0f)
{
    timestamp_prev = _micros();
}


float LowPassFilter::operator() (float x)
{
    unsigned long timestamp = _micros();
    float dt = (timestamp - timestamp_prev)*1e-6f;

    if (dt < 0.0f || dt > 0.5f)
        dt = 1e-3f;

    float alpha = Tf/(Tf + dt);
    float y = alpha*y_prev + (1.0f - alpha)*x;

    y_prev = y;
    timestamp_prev = timestamp;
    return y;
}
```

### 4）pid

> #### 头文件类定义

```c++
/**
 *  PID controller class
 */
class PIDController
{
public:
    /**
     *  
     * @param P - Proportional gain 
     * @param I - Integral gain
     * @param D - Derivative gain 
     * @param ramp - Maximum speed of change of the output value
     * @param limit - Maximum output value
     */
    PIDController(float P, float I, float D, float ramp, float limit);
    ~PIDController() = default;

    float operator() (float error);

    float P; //!< Proportional gain 
    float I; //!< Integral gain 
    float D; //!< Derivative gain 
    float output_ramp; //!< Maximum speed of change of the output value
    float limit; //!< Maximum output value

protected:
    float integral_prev; //!< last integral component value
    float error_prev; //!< last tracking error value
    float timestamp_prev; //!< Last execution timestamp
    float output_prev;  //!< last pid output value
};
```

> ####  具体函数

```c++
PIDController::PIDController(float P, float I, float D, float ramp, float limit)
    : P(P)
    , I(I)
    , D(D)
    , output_ramp(ramp)    // 输出导数限制（斜坡） [volts/second]
    , limit(limit)         // 输出数值（电压）限制 [volts]
    , integral_prev(0.0)
    , error_prev(0.0)
    , output_prev(0.0)
{
    timestamp_prev = _micros()*1e-6;
}

// PID控制器功能
float PIDController::operator() (float error){
    // 计算上一次时间
    unsigned long timestamp_now = _micros();
    float Ts = (timestamp_now - timestamp_prev) * 1e-6;
    // 异常情况修复（时间溢出）
    if(Ts <= 0 || Ts > 0.5) Ts = 1e-3; 

    // u(s) = (P + I/s + Ds)e(s)
    // 实现离散
    // 比例部分 
    // u_p  = P *e(k)
    float proportional = P * error;
    // 积分部分双线性变换
    // u_ik = u_ik_1  + I*Ts/2*(ek + ek_1)
    float integral = integral_prev + I*Ts*0.5*(error + error_prev);
    // 抗饱和-限制输出电压
    integral = _constrain(integral, -limit, limit);
    // 离散导数
    // u_dk = D(ek - ek_1)/Ts
    float derivative = D*(error - error_prev)/Ts;

    // 所有部分总和
    float output = proportional + integral + derivative;
    // 限制输出变量
    output = _constrain(output, -limit, limit);

    // 斜坡处理，限制输出加速度
    float output_rate = (output - output_prev)/Ts;
    if (output_rate > output_ramp)
        output = output_prev + output_ramp*Ts;
    else if (output_rate < -output_ramp)
        output = output_prev - output_ramp*Ts;

    // 保存本次变量
    integral_prev = integral;
    output_prev = output;
    error_prev = error;
    timestamp_prev = timestamp_now;
    
    return output;
}
```

### 5）time_utils

主要包含延时函数（ _delay）和获取时间戳函数（ _micros）



### 6) base_classes

#### a. 无刷直流电机驱动

```C++
class BLDCDriver{
    public:       
        /** 初始化硬件 */
        virtual int init();
        /** 开启硬件 */
        virtual void enable();
        /** 关闭硬件 */
        virtual void disable();
        long pwm_frequency; //PWM频率（hz）
        float voltage_power_supply; //电源电压
        float voltage_limit; //限制电机电压    
        /** 
         * 设置相电压 （ABC相）
         * 
         * @param Ua - phase A voltage
         * @param Ub - phase B voltage
         * @param Uc - phase C voltage
        */
        virtual void setPwm(float Ua, float Ub, float Uc);
};
```

#### b. 步进电机驱动

```C++
class StepperDriver{
    public:
        
        /** 初始化硬件 */
        virtual int init();
        /** 开启 */
        virtual void enable();
        /** 关闭 */
        virtual void disable();

        long pwm_frequency; //pwm频率（hz）
        float voltage_power_supply; //电源电压
        float voltage_limit; //限制电机电压
            
        /** 
         * 设置相电压 （AB相）
         * 
         * @param Ua phase A voltage
         * @param Ub phase B voltage
        */
        virtual void setPwm(float Ua, float Ub);
};
```

#### c. 传感器（编码器）

```c++
/**
 *  方向枚举
 */
enum Direction{
    CW      = 1,  //clockwise 顺时针
    CCW     = -1, //counter clockwise 逆时针
    UNKNOWN = 0   //无知或无效状态
};

/**
 *  上拉配置枚举
 */
enum Pullup{
    INTERN, //使用内部上拉
    EXTERN  //使用外部上拉
};

/**
 *  传感器抽象类定义
 *  每个传感器都需要实现这些功能
 */
class Sensor{
    public:
        /** 获取当前角度 (rad) */
        virtual float getAngle();
        /** 获取当前角速度 (rad/s)*/
        virtual float getVelocity();
        /**
         *  将当前角度设置为零角度
         *  返回角度差（rad）
         */
        virtual float initRelativeZero();
        /**
         * 将绝对零角设置为零角
         * 返回角度差（rad）
         */
        virtual float initAbsoluteZero();

        // 如果自然方向为ccw逆时针，则将方向翻转为cw顺时针
        int natural_direction = Direction::CW;

        /** 
         * 如果没有绝对0度值，则返回0
         * 0 - 无标志增量式编码器
         * 1 - 带标志和磁传感器的编码器
         */
        virtual int hasAbsoluteZero();
        /** 
         * 如果要搜索绝对0度值，则返回0
         * 0 - 磁传感器 (找到有标志的编码器)
         * 1 - 带标志的编码器 (尚未找到标志)
         */
        virtual int needsAbsoluteZeroSearch();
};
```

#### d. FOC相关

##### 枚举、变量定义及函数定义说明

```C++

/**
 *  运动控制类型枚举
 */
enum ControlType{
  voltage,			//电压转矩闭环控制
  velocity,			//速度闭环控制
  angle,			//位置/角度闭环控制
  velocity_openloop, //速度开环控制
  angle_openloop	//角度开环控制
};

/**
 *  FOC调制类型
 */
enum FOCModulationType{
  SinePWM, 			//正弦PWM调制
  SpaceVectorPWM, 	//空间矢量脉宽调制（SVPWM）
  Trapezoid_120,  	//梯形轨迹120
  Trapezoid_150	  	//梯形轨迹150
};

/**
 通用电机类
*/
class FOCMotor
{
  public:
    /**
     * 默认构造函数，将所有变量设置为默认值
     */
    FOCMotor();

    /**  电机硬件初始化 */
  	virtual void init()=0;
    /** 电机禁用 */
  	virtual void disable()=0;
    /** 电机使能 */
    virtual void enable()=0;

    /**
     * 连接电机和传感器 
     * @param  Sensor sensor类，用于FOC算法读取电机角度和速度
     */
    void linkSensor(Sensor* sensor);

    
    /**
     * 初始化FOC算法函数
     * 调整传感器和电机的零位 
     * 如果设置了zero_electric_offset 参数，则跳过对齐过程
     * 
     * @param zero_electric_offset 相对电机0位置的传感器的绝对位置偏移值.
     * @param sensor_direction  传感器自然转向，默认为CW顺时针
     *
     */  
    virtual int initFOC( float zero_electric_offset = NOT_SET , Direction sensor_direction = Direction::CW)=0; 
    
    /**
     * 函数实时运行FOC算法
     * 计算电机角度并设置适合的电压 
     * 相位PWM信号
     * 运行周期越快越好
     */ 
    virtual void loopFOC()=0;
    
    /**
     * 执行由无刷直流电机控制器参数设置的控制回路的函数.
     * @param target  基于电机驱动的目标电压、角度或速度
     *  如果这个参数未设置，电机将使用其变量中设置的目标运动目标
     *  这个函数不需要在每次循环执行时运行，取决于如何使用
     */
    virtual void move(float target = NOT_SET)=0;

    // 状态计算方法 
    /** 轴角度计算 [rad] */
    float shaftAngle();
    /** 
     * 轴角计算函数 [rad/s]
     * 实现了低通滤波
     */
    float shaftVelocity();

    // 状态变量
    float target; //当前目标值
  	float shaft_angle;//当前电机角度
  	float shaft_velocity;//当前电机速度
    float shaft_velocity_sp;//当前目标速度
    float shaft_angle_sp;//当前目标角度
 	float voltage_q;//当前 u_q 设定
    float voltage_d;//当前 u_d 设定
    float voltage_sensor_align;//传感器和电机对齐电压参数
    float velocity_index_search;//标志搜索的目标速度
    int pole_pairs;	//电机极对数

    //限制变量
    float voltage_limit; //全局电压限制
    float velocity_limit; //全局速度限制

    float zero_electric_angle;//绝对零角限制（如果有的话）

    // 初始化构造
    ControlType controller; //确定要使用的控制回路参数
    FOCModulationType foc_modulation;//foc算法
    PIDController PID_velocity{DEF_PID_VEL_P,DEF_PID_VEL_I,DEF_PID_VEL_D,DEF_PID_VEL_U_RAMP,DEF_POWER_SUPPLY};//速度环PI配置参数
    PIDController P_angle{DEF_P_ANGLE_P,0,0,1e10,DEF_VEL_LIM};//位置环P配置参数
    LowPassFilter LPF_velocity{DEF_VEL_FILTER_Tf};//速度低通滤波器配置参数

    /**
     * 函数为BLDCMotor类提供串口和启用监控功能
     * @param serial 
     */
    void useMonitoring(Print &serial);

    /**
     * 与串口绘图仪一起使用的实用功能，用于监控电机变量，这会大大降低执行速度！
     */
    void monitor();

 /**
*功能设定电机的配置参数，控制回路的目标值
*并将其输出到监控端口（如果可用）：
*-配置PID控制器常数
*-更改运动控制回路
*-监控电机变量
*-设定目标值
*-检查所有配置值
*
*要检查配置值，只需输入命令字母。
*例如：
*-读取速度PI控制器P增益运行：P
*-将速度PI控制器P增益设置为1.2运行：P1.2
*
*要更改目标值，只需在终端中输入一个数字：
*例如：
*-要将目标值更改为-0.1453，请输入：-0.1453
*-要获取当前目标值，请输入：V3
*
*命令列表：
*-P：速度PI控制器P增益
*I：速度PI控制器I增益
*L：速度PI控制器电压限制
*R：速度PI控制器电压斜坡
*-F：速度低通滤波器时间常数
*K：角度P控制器P增益
*N：P角控制器速度限制
*C：控制回路
*-0：电压
*-1:速度
*-2：角度
*-V:获取电机变量
*-0：当前设定电压
*-1：流速
*-2：电流角
*-3：当前目标值
*@param command包含用户命令的字符串
*对于错误返回0，对于已执行的命令返回1
*/
    int command(String command);
    
    /** 
      * 连接传感器:
      * - 编码器 
      * - 磁编码器
      * - 霍尔传感器
    */
    Sensor* sensor; 
    // 监控功能
    Print* monitor_port; //串口终端变量（如果有）
};

```

##### 具体的函数定义：

```C++

/**
 * Default constructor - setting all variabels to default values
 */
FOCMotor::FOCMotor()
{
  // maximum angular velocity to be used for positioning 
  velocity_limit = DEF_VEL_LIM;
  // maximum voltage to be set to the motor
  voltage_limit = DEF_POWER_SUPPLY;

  // index search velocity
  velocity_index_search = DEF_INDEX_SEARCH_TARGET_VELOCITY;
  // sensor and motor align voltage
  voltage_sensor_align = DEF_VOLTAGE_SENSOR_ALIGN;

  // electric angle of comthe zero angle
  zero_electric_angle = 0;

  // default modulation is SinePWM
  foc_modulation = FOCModulationType::SinePWM;

  // default target value
  target = 0;
  voltage_d = 0;
  voltage_q = 0;
  
  //monitor_port 
  monitor_port = nullptr;
  //sensor 
  sensor = nullptr;
}


/**
	Sensor communication methods
*/
void FOCMotor::linkSensor(Sensor* _sensor) {
  sensor = _sensor;
}
// shaft angle calculation
float FOCMotor::shaftAngle() {
  // if no sensor linked return 0
  if(!sensor) return shaft_angle;
  return sensor->getAngle();
}
// shaft velocity calculation
float FOCMotor::shaftVelocity() {
  // if no sensor linked return 0
  if(!sensor) return 0;
  return LPF_velocity(sensor->getVelocity());
}

/**
 *  Monitoring functions
 */
// function implementing the monitor_port setter
void FOCMotor::useMonitoring(Print &print){
  monitor_port = &print; //operate on the address of print
  if(monitor_port ) monitor_port->println("MOT: Monitor enabled!");
}
// utility function intended to be used with serial plotter to monitor motor variables
// significantly slowing the execution down!!!!
void FOCMotor::monitor() {
  if(!monitor_port) return;
  switch (controller) {
    case ControlType::velocity_openloop:
    case ControlType::velocity:
      monitor_port->print(voltage_q);
      monitor_port->print("\t");
      monitor_port->print(shaft_velocity_sp);
      monitor_port->print("\t");
      monitor_port->println(shaft_velocity);
      break;
    case ControlType::angle_openloop:
    case ControlType::angle:
      monitor_port->print(voltage_q);
      monitor_port->print("\t");
      monitor_port->print(shaft_angle_sp);
      monitor_port->print("\t");
      monitor_port->println(shaft_angle);
      break;
    case ControlType::voltage:
      monitor_port->print(voltage_q);
      monitor_port->print("\t");
      monitor_port->print(shaft_angle);
      monitor_port->print("\t");
      monitor_port->println(shaft_velocity);
      break;
  }
}

int FOCMotor::command(String user_command) {
  // error flag
  int errorFlag = 1;
  // if empty string
  if(user_command.length() < 1) return errorFlag;

  // parse command letter
  char cmd = user_command.charAt(0);
  // check if get command
  char GET = user_command.charAt(1) == '\n';
  // parse command values
  float value = user_command.substring(1).toFloat();

  // a bit of optimisation of variable memory for Arduino UNO (atmega328)
  switch(cmd){
    case 'P':      // velocity P gain change
    case 'I':      // velocity I gain change
    case 'D':      // velocity D gain change
    case 'R':      // velocity voltage ramp change
      if(monitor_port) monitor_port->print(" PID velocity| ");
      break;
    case 'F':      // velocity Tf low pass filter change
      if(monitor_port) monitor_port->print(" LPF velocity| ");
      break;
    case 'K':      // angle loop gain P change
      if(monitor_port) monitor_port->print(" P angle| ");
      break;
    case 'L':      // velocity voltage limit change
    case 'N':      // angle loop gain velocity_limit change
      if(monitor_port) monitor_port->print(" Limits| ");
      break;
  }

  // apply the the command
  switch(cmd){
    case 'P':      // velocity P gain change
      if(monitor_port) monitor_port->print("P: ");
      if(!GET) PID_velocity.P = value;
      if(monitor_port) monitor_port->println(PID_velocity.P);
      break;
    case 'I':      // velocity I gain change
      if(monitor_port) monitor_port->print("I: ");
      if(!GET) PID_velocity.I = value;
      if(monitor_port) monitor_port->println(PID_velocity.I);
      break;
    case 'D':      // velocity D gain change
      if(monitor_port) monitor_port->print("D: ");
      if(!GET) PID_velocity.D = value;
      if(monitor_port) monitor_port->println(PID_velocity.D);
      break;
    case 'R':      // velocity voltage ramp change
      if(monitor_port) monitor_port->print("volt_ramp: ");
      if(!GET) PID_velocity.output_ramp = value;
      if(monitor_port) monitor_port->println(PID_velocity.output_ramp);
      break;
    case 'L':      // velocity voltage limit change
      if(monitor_port) monitor_port->print("volt_limit: ");
      if(!GET) {
        voltage_limit = value;
        PID_velocity.limit = value;
      }
      if(monitor_port) monitor_port->println(voltage_limit);
      break;
    case 'F':      // velocity Tf low pass filter change
      if(monitor_port) monitor_port->print("Tf: ");
      if(!GET) LPF_velocity.Tf = value;
      if(monitor_port) monitor_port->println(LPF_velocity.Tf);
      break;
    case 'K':      // angle loop gain P change
      if(monitor_port) monitor_port->print(" P: ");
      if(!GET) P_angle.P = value;
      if(monitor_port) monitor_port->println(P_angle.P);
      break;
    case 'N':      // angle loop gain velocity_limit change
      if(monitor_port) monitor_port->print("vel_limit: ");
      if(!GET){
        velocity_limit = value;
        P_angle.limit = value;
      }
      if(monitor_port) monitor_port->println(velocity_limit);
      break;
    case 'C':
      // change control type
      if(monitor_port) monitor_port->print("Control: ");
      
      if(GET){ // if get command
        switch(controller){
          case ControlType::voltage:
            if(monitor_port) monitor_port->println("voltage");
            break;
          case ControlType::velocity:
            if(monitor_port) monitor_port->println("velocity");
            break;
          case ControlType::angle:
            if(monitor_port) monitor_port->println("angle");
            break;
          default:
            if(monitor_port) monitor_port->println("open loop");
        }
      }else{ // if set command
        switch((int)value){
          case 0:
            if(monitor_port) monitor_port->println("voltage");
            controller = ControlType::voltage;
            break;
          case 1:
            if(monitor_port) monitor_port->println("velocity");
            controller = ControlType::velocity;
            break;
          case 2:
            if(monitor_port) monitor_port->println("angle");
            controller = ControlType::angle;
            break;
          default: // not valid command
            errorFlag = 0;
        }
      }
      break;
    case 'V':     // get current values of the state variables
        switch((int)value){
          case 0: // get voltage
            if(monitor_port) monitor_port->print("Uq: ");
            if(monitor_port) monitor_port->println(voltage_q);
            break;
          case 1: // get velocity
            if(monitor_port) monitor_port->print("Velocity: ");
            if(monitor_port) monitor_port->println(shaft_velocity);
            break;
          case 2: // get angle
            if(monitor_port) monitor_port->print("Angle: ");
            if(monitor_port) monitor_port->println(shaft_angle);
            break;
          case 3: // get angle
            if(monitor_port) monitor_port->print("Target: ");
            if(monitor_port) monitor_port->println(target);
            break;
          default: // not valid command
            errorFlag = 0;
        }
      break;
    default:  // target change
      if(monitor_port) monitor_port->print("Target : ");
      target = user_command.toFloat();
      if(monitor_port) monitor_port->println(target);
  }
  // return 0 if error and 1 if ok
  return errorFlag;
}
```





## 2.drivers	驱动文件

### 1）hardware_specific

特定硬件配置，为了适配不同mcu而使用不同的配置

> atmega328_mcu
>
> atmega2560_mcu
>
> esp32_mcu
>
> stm32_mcu
>
> teensy_mcu

### 2) BLDCDriver3PWM

> 头文件：

```c++
/**
 3 pwm bldc 驱动器类
*/
class BLDCDriver3PWM: public BLDCDriver
{
  public:
    /**
      BLDCDriver class constructor
      @param phA A phase pwm pin
      @param phB B phase pwm pin
      @param phC C phase pwm pin
      @param en enable pin (optional input)
    */
    BLDCDriver3PWM(int phA,int phB,int phC, int en = NOT_SET);
    
    /**  Motor hardware init function */
  	int init() override;
    /** Motor disable function */
  	void disable() override;
    /** Motor enable function */
    void enable() override;

    // 变量
  	int pwmA; //!< phase A pwm pin number
  	int pwmB; //!< phase B pwm pin number
  	int pwmC; //!< phase C pwm pin number
    int enable_pin; //!< enable pin number

    /** 
     * 设定硬件相电压
     * 
     * @param Ua - phase A voltage
     * @param Ub - phase B voltage
     * @param Uc - phase C voltage
    */
    void setPwm(float Ua, float Ub, float Uc) override;

  private:
};

```

> 初始化：

```c++
BLDCDriver3PWM::BLDCDriver3PWM(int phA, int phB, int phC, int en){
  // 引脚初始化
  pwmA = phA;
  pwmB = phB;
  pwmC = phC;
  // 使能引脚
  enable_pin = en;
  // 默认电源值
  voltage_power_supply = DEF_POWER_SUPPLY;
  voltage_limit = NOT_SET;
}

// 使能电机驱动
void  BLDCDriver3PWM::enable(){
    // enable_pin the driver - if enable_pin pin available
    if ( enable_pin != NOT_SET ) digitalWrite(enable_pin, HIGH);
    // set zero to PWM
    setPwm(0,0,0);
}

// 禁用电机驱动
void BLDCDriver3PWM::disable()
{
  // set zero to PWM
  setPwm(0, 0, 0);
  // disable the driver - if enable_pin pin available
  if ( enable_pin != NOT_SET ) digitalWrite(enable_pin, LOW);

}

// 初始化硬件针脚   
int BLDCDriver3PWM::init() {
  // a bit of separation
  _delay(1000);

  // PWM pins
  pinMode(pwmA, OUTPUT);
  pinMode(pwmB, OUTPUT);
  pinMode(pwmC, OUTPUT);
  if(enable_pin != NOT_SET) pinMode(enable_pin, OUTPUT);

  // 电压限制配置的健全性检查
  if(voltage_limit == NOT_SET || voltage_limit > voltage_power_supply) voltage_limit =  voltage_power_supply;

  // 将pwm频率设置为引脚
  // 硬件特定功能-取决于驱动程序和mcu
  _configure3PWM(pwm_frequency, pwmA, pwmB, pwmC);
  return 0;
}
```

> 设定pwm引脚相电压

```c++
void BLDCDriver3PWM::setPwm(float Ua, float Ub, float Uc) {  
  // 限制驱动器电压
  Ua = _constrain(Ua, 0.0, voltage_limit);
  Ub = _constrain(Ub, 0.0, voltage_limit);
  Uc = _constrain(Uc, 0.0, voltage_limit);    
  // 计算占空比
  // 限制在 [0,1]
  float dc_a = _constrain(Ua / voltage_power_supply, 0.0 , 1.0 );
  float dc_b = _constrain(Ub / voltage_power_supply, 0.0 , 1.0 );
  float dc_c = _constrain(Uc / voltage_power_supply, 0.0 , 1.0 );
  // 硬件特定写入
  // 硬件特定功能-取决于驱动程序和mcu
  _writeDutyCycle3PWM(dc_a, dc_b, dc_c, pwmA, pwmB, pwmC);
}
```

### 3）BLDCDriver6PWM

跟3pwm类似配置

> 头文件

```c++
/**
 6 pwm bldc driver class
*/
class BLDCDriver6PWM: public BLDCDriver
{
  public:
    /**
      BLDCDriver class constructor
      @param phA_h A phase pwm pin
      @param phA_l A phase pwm pin
      @param phB_h B phase pwm pin
      @param phB_l A phase pwm pin
      @param phC_h C phase pwm pin
      @param phC_l A phase pwm pin
      @param en enable pin (optional input)
    */
    BLDCDriver6PWM(int phA_h,int phA_l,int phB_h,int phB_l,int phC_h,int phC_l, int en = NOT_SET);
    
    /**  Motor hardware init function */
  	int init() override;
    /** Motor disable function */
  	void disable() override;
    /** Motor enable function */
    void enable() override;

    // hardware variables
  	int pwmA_h,pwmA_l; //!< phase A pwm pin number
  	int pwmB_h,pwmB_l; //!< phase B pwm pin number
  	int pwmC_h,pwmC_l; //!< phase C pwm pin number
    int enable_pin; //!< enable pin number

    float dead_zone; //!< a percentage of dead-time(zone) (both high and low side in low) for each pwm cycle [0,1]

    /** 
     * Set phase voltages to the harware 
     * 
     * @param Ua - phase A voltage
     * @param Ub - phase B voltage
     * @param Uc - phase C voltage
    */
    void setPwm(float Ua, float Ub, float Uc) override;

  private:
        
};

```

> 函数构建：

```c++
BLDCDriver6PWM::BLDCDriver6PWM(int phA_h,int phA_l,int phB_h,int phB_l,int phC_h,int phC_l, int en){
  // Pin initialization
  pwmA_h = phA_h;
  pwmB_h = phB_h;
  pwmC_h = phC_h;
  pwmA_l = phA_l;
  pwmB_l = phB_l;
  pwmC_l = phC_l;

  // enable_pin pin
  enable_pin = en;

  // default power-supply value
  voltage_power_supply = DEF_POWER_SUPPLY;
  voltage_limit = NOT_SET;
 
  // dead zone initial - 2%
  dead_zone = 0.02;

}

// enable motor driver
void  BLDCDriver6PWM::enable(){
    // enable_pin the driver - if enable_pin pin available
    if ( enable_pin != NOT_SET ) digitalWrite(enable_pin, HIGH);
    // set zero to PWM
    setPwm(0, 0, 0);
}

// disable motor driver
void BLDCDriver6PWM::disable()
{
  // set zero to PWM
  setPwm(0, 0, 0);
  // disable the driver - if enable_pin pin available
  if ( enable_pin != NOT_SET ) digitalWrite(enable_pin, LOW);

}

// init hardware pins   
int BLDCDriver6PWM::init() {
  // a bit of separation
  _delay(1000);

  // PWM pins
  pinMode(pwmA_l, OUTPUT);
  pinMode(pwmB_h, OUTPUT);
  pinMode(pwmC_h, OUTPUT);
  pinMode(pwmA_l, OUTPUT);
  pinMode(pwmB_l, OUTPUT);
  pinMode(pwmC_l, OUTPUT);
  if(enable_pin != NOT_SET) pinMode(enable_pin, OUTPUT);


  // sanity check for the voltage limit configuration
  if(voltage_limit == NOT_SET || voltage_limit > voltage_power_supply) voltage_limit =  voltage_power_supply;

  // configure 6pwm 
  // hardware specific function - depending on driver and mcu
  return _configure6PWM(pwm_frequency, dead_zone, pwmA_h,pwmA_l, pwmB_h,pwmB_l, pwmC_h,pwmC_l);
}


// Set voltage to the pwm pin
void BLDCDriver6PWM::setPwm(float Ua, float Ub, float Uc) {  
  // limit the voltage in driver
  Ua = _constrain(Ua, 0, voltage_limit);
  Ub = _constrain(Ub, 0, voltage_limit);
  Uc = _constrain(Uc, 0, voltage_limit);    
  // calculate duty cycle
  // limited in [0,1]
  float dc_a = _constrain(Ua / voltage_power_supply, 0.0 , 1.0 );
  float dc_b = _constrain(Ub / voltage_power_supply, 0.0 , 1.0 );
  float dc_c = _constrain(Uc / voltage_power_supply, 0.0 , 1.0 );
  // hardware specific writing
  // hardware specific function - depending on driver and mcu
  _writeDutyCycle6PWM(dc_a, dc_b, dc_c, dead_zone, pwmA_h,pwmA_l, pwmB_h,pwmB_l, pwmC_h,pwmC_l);
}
```

### 4）StepperDriver2PWM

### 5）StepperDriver4PWM

### 6）hardware_api

#### a. 初始化

```c++
/** 
 * 配置PWM频率、分辨率和校准
 * - 步进驱动器-2PWM设置
 * - hardware specific
 * 
 * @param pwm_frequency - 频率（赫兹）-如果适用
 * @param pinA pinA bldc驱动
 * @param pinB pinB bldc驱动
 */
void _configure2PWM(long pwm_frequency, const int pinA, const int pinB);

/** 
 * 配置PWM频率、分辨率和校准
 * - BLDC驱动器-3PWM设置
 * - hardware specific
 * 
 * @param pwm_frequency - frequency in hertz - if applicable
 * @param pinA pinA bldc driver
 * @param pinB pinB bldc driver
 * @param pinC pinC bldc driver
 */
void _configure3PWM(long pwm_frequency, const int pinA, const int pinB, const int pinC);

/** 
 * 配置PWM频率、分辨率和校准
 * - 步进驱动器-4PWM设置
 * - hardware specific
 * 
 * @param pwm_frequency - frequency in hertz - if applicable
 * @param pin1A pin1A stepper driver
 * @param pin1B pin1B stepper driver
 * @param pin2A pin2A stepper driver
 * @param pin2B pin2B stepper driver
 */
void _configure4PWM(long pwm_frequency, const int pin1A, const int pin1B, const int pin2A, const int pin2B);

/** 
 * 配置PWM频率、分辨率和校准
 * - BLDC驱动器-6PWM设置
 * - hardware specific
 * 
 * @param pwm_frequency - frequency in hertz - if applicable
 * @param dead_zone  duty cycle protection zone [0, 1] - both low and high side low - if applicable
 * @param pinA_h pinA high-side bldc driver 
 * @param pinA_l pinA low-side bldc driver 
 * @param pinB_h pinA high-side bldc driver 
 * @param pinB_l pinA low-side bldc driver 
 * @param pinC_h pinA high-side bldc driver 
 * @param pinC_l pinA low-side bldc driver 
 * 
 * @return 0 if config good, -1 if failed
 */
int _configure6PWM(long pwm_frequency, float dead_zone, const int pinA_h, const int pinA_l,  const int pinB_h, const int pinB_l, const int pinC_h, const int pinC_l);

```

#### b. 设置占空比

```c++
/** 
 * 设置pwm引脚的占空比功能 (ex. analogWrite())
 * - 步进驱动器-2PWM设置
 * - hardware specific
 * 
 * @param dc_a  duty cycle phase A [0, 1]
 * @param dc_b  duty cycle phase B [0, 1]
 * @param pinA  phase A hardware pin number
 * @param pinB  phase B hardware pin number
 */ 
void _writeDutyCycle2PWM(float dc_a,  float dc_b, int pinA, int pinB);

/** 
 * 设置pwm引脚的占空比功能 (ex. analogWrite())
 * - BLDC 驱动器 - 3PWM设置
 * - hardware specific
 * 
 * @param dc_a  duty cycle phase A [0, 1]
 * @param dc_b  duty cycle phase B [0, 1]
 * @param dc_c  duty cycle phase C [0, 1]
 * @param pinA  phase A hardware pin number
 * @param pinB  phase B hardware pin number
 * @param pinC  phase C hardware pin number
 */ 
void _writeDutyCycle3PWM(float dc_a,  float dc_b, float dc_c, int pinA, int pinB, int pinC);

/** 
 * 设置pwm引脚的占空比功能 (ex. analogWrite())
 * - 步进驱动器 - 4PWM setting
 * - hardware specific
 * 
 * @param dc_1a  duty cycle phase 1A [0, 1]
 * @param dc_1b  duty cycle phase 1B [0, 1]
 * @param dc_2a  duty cycle phase 2A [0, 1]
 * @param dc_2b  duty cycle phase 2B [0, 1]
 * @param pin1A  phase 1A hardware pin number
 * @param pin1B  phase 1B hardware pin number
 * @param pin2A  phase 2A hardware pin number
 * @param pin2B  phase 2B hardware pin number
 */ 
void _writeDutyCycle4PWM(float dc_1a,  float dc_1b, float dc_2a, float dc_2b, int pin1A, int pin1B, int pin2A, int pin2B);


/** 
 * 设置pwm引脚的占空比功能 (ex. analogWrite())
 * - BLDC 驱动器 - 6PWM setting
 * - hardware specific
 * 
 * @param dc_a  duty cycle phase A [0, 1]
 * @param dc_b  duty cycle phase B [0, 1]
 * @param dc_c  duty cycle phase C [0, 1]
 * @param dead_zone  duty cycle protection zone [0, 1] - both low and high side low
 * @param pinA_h  phase A high-side hardware pin number
 * @param pinA_l  phase A low-side hardware pin number
 * @param pinB_h  phase B high-side hardware pin number
 * @param pinB_l  phase B low-side hardware pin number
 * @param pinC_h  phase C high-side hardware pin number
 * @param pinC_l  phase C low-side hardware pin number
 * 
 */ 
void _writeDutyCycle6PWM(float dc_a,  float dc_b, float dc_c, float dead_zone, int pinA_h, int pinA_l, int pinB_h, int pinB_l, int pinC_h, int pinC_l);

```

## 3.sensors    