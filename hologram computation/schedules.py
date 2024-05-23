# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 15:50:22 2022

@author: asupikov

Almost direct port of hypernerf schedule system
"""

def create_schedule(config):
    sched = None
    if config["otype"]=="none":
        sched = NoneSchedule()
    elif config["otype"]=="constant":
        sched = ConstantSchedule(config["val"])
    elif config["otype"] == "linear":
        sched = LinearSchedule(config["init_val"], config["final_val"], config["num_steps"])
    elif config["otype"] == "exponential":
        sched = ExponentialSchedule(config["init_val"], config["final_val"], config["num_steps"])
    if("delay_steps" in config):
        return DelayedSchedule(sched, config["delay_steps"])    
    return sched
   
    
class Schedule:
  """An interface for generic schedules.."""

  def get(self, step):
    """Get the value for the given step."""
    raise NotImplementedError

  def __call__(self, step):
    return self.get(step)


class NoneSchedule(Schedule):
  """Always returns None. Useful for disable schedules."""
  def __init__(self):
    super().__init__()

  def get(self, step):
    return None


class ConstantSchedule(Schedule):
  """Linearly scaled scheduler."""

  def __init__(self, value):
    super().__init__()
    self.value = value

  def get(self, step):
    """Get the value for the given step."""
    if self.value is None:
      return None
    return self.value


class LinearSchedule(Schedule):
  """Linearly scaled scheduler."""

  def __init__(self, initial_value, final_value, num_steps):
    super().__init__()
    self.initial_value = initial_value
    self.final_value = final_value
    self.num_steps = num_steps
    
  def get(self, step):
    """Get the value for the given step."""
    
    if self.num_steps == 0:
      return self.final_value
    alpha = min(step / self.num_steps, 1.0)
    return (1.0 - alpha) * self.initial_value + alpha * self.final_value


class ExponentialSchedule(Schedule):
  """Exponentially decaying scheduler."""

  def __init__(self, initial_value, final_value, num_steps, eps=1e-10):
    super().__init__()
    if initial_value <= final_value:
      raise ValueError('Final value must be less than initial value.')

    self.initial_value = initial_value
    self.final_value = final_value
    self.num_steps = num_steps
    self.eps = eps

  def get(self, step):
    """Get the value for the given step."""
    if step >= self.num_steps:
      return self.final_value

    final_value = max(self.final_value, self.eps)
    base = final_value / self.initial_value
    exponent = step / (self.num_steps - 1)
    if step >= self.num_steps:
      return self.final_value
    return self.initial_value * base**exponent

class DelayedSchedule(Schedule) :
    def __init__(self, schedule, delay_steps):
      super().__init__()
      self.delay = delay_steps
      self.schedule = schedule
      
    def get(self, step) :
        if step >= self.delay :
            return self.schedule(step-self.delay)
        return self.schedule(0)
        
