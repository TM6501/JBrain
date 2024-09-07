# This file is meant to be used from within C++ and is thus written
# to be as easy-to-use as possible from the clunky C++ Python API.
# As much as possible is hard-coded.

import torch
import math
import numpy as np
import gymnasium as gym

from stable_baselines3 import PPO
from stable_baselines3.common.env_util import make_vec_env, DummyVecEnv
from stable_baselines3.common.callbacks import EvalCallback, StopTrainingOnRewardThreshold

class SageRunner:
	def __init__(self):
		self.useIntAction = True	
	
	# We want a simple 1-environment agent:
	def getEnv(self, renderMode=None):
		return DummyVecEnv([lambda: gym.make("CartPole-v1", render_mode=renderMode)])

	# Training is independent and won't be called from C++:
	def train(self):
		vec_env = getEnv()
		eval_env = gym.make("CartPole-v1")
	
		callback_on_best = StopTrainingOnRewardThreshold(reward_threshold=500, verbose=1)
		eval_callback = EvalCallback(eval_env, callback_on_new_best=callback_on_best, verbose=1)	

		model = PPO("MlpPolicy", vec_env, verbose=1)

		model.learn(total_timesteps=100000, callback=eval_callback)
		model.save("ppo_sage")

	def load(self):

		self.model = PPO.load('ppo_sage')
		self.env = self.getEnv(renderMode='rgb_array')

	def render(self):		
		self.env.render()

	def reset(self):
		# Gather observations and actions:
		obs = self.env.reset()
		actions, states = self.model.predict(obs)
	
		# We only want the first observation and first action:
		actions = self.getList(actions[0])
		obs = self.getList(obs[0])
		
		return obs, actions

	# We often need to make sure our return value is a list. Do that here:
	def getList(self, inputVal):
		if isinstance(inputVal, np.ndarray):
			return list(inputVal)
		else:
			return list([inputVal])

	# Reminder, as a sage it wants to return an observation and what the sage-model
	# would do under those circumstances.
	def step(self, actualAction):
		if self.useIntAction:
			actualAction = int(actualAction)
		
		inAction = (actualAction, )
		obs, rewards, dones, info = self.env.step(inAction)
		actions, states = self.model.predict(obs)
	
		obs = self.getList(obs[0])
		actions = self.getList(actions[0])
		rewards = rewards[0]
		dones = dones[0]
				
		# Vec-env returns need to be converted to what the singular
		# agent needs:
		return (obs, float(rewards), bool(dones), actions)

	def test(self):
		vec_env = self.getEnv(renderMode='human')
		model = PPO.load("ppo_sage")

		obs = vec_env.reset()
		maxSteps = 2000
		# With a simple environment like this, what we want is in subscript
		# 0 of each of action, rewards, dones, and info.
		for i in range(maxSteps):
			action, _states = model.predict(obs)
			obs, rewards, dones, info = vec_env.step(action)
			vec_env.render()
			if dones[0]:
				print(f"Step {i}: Done. Resetting the environment.")

def getInstance():
	return SageRunner()

if __name__ == "__main__":
	# train()
	# test()
	print("In Python main.", flush=True)